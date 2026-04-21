import csv

import pytest

from src import map_variants as mv


def _write_tsv(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["variant_urn", "raw_hgvs_nt", "raw_hgvs_pro", "target_sequence"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_map_variants_default_preserves_input_order(tmp_path, monkeypatch):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"

    rows = [
        {"variant_urn": "v0", "raw_hgvs_nt": "", "raw_hgvs_pro": "p.Ala1Val", "target_sequence": "SEQ_A"},
        {"variant_urn": "v1", "raw_hgvs_nt": "", "raw_hgvs_pro": "p.Ala2Val", "target_sequence": "SEQ_B"},
        {"variant_urn": "v2", "raw_hgvs_nt": "", "raw_hgvs_pro": "p.Ala3Val", "target_sequence": "SEQ_A"},
        {"variant_urn": "v3", "raw_hgvs_nt": "", "raw_hgvs_pro": "p.Ala4Val", "target_sequence": "SEQ_B"},
    ]
    _write_tsv(input_path, rows)

    async def fake_pipeline(group_name, target_seq, row_entries, dcd):
        per_row = []
        for orig_idx, _, raw_pro, _ in row_entries:
            per_row.append((orig_idx, f"NC_000001.11:g.{orig_idx + 100}A>G", None))
        return per_row, "NM_000001.1"

    async def fake_clingen_batch(hgvs_strings, max_concurrency=5):
        return {h: {"hgvs": h, "id": f"CA{idx}"} for idx, h in enumerate(hgvs_strings, start=1)}

    monkeypatch.setattr(mv, "_try_import_dcd_mapping", lambda: object())
    monkeypatch.setattr(mv, "_run_dcd_mapping_pipeline", fake_pipeline)
    monkeypatch.setattr(mv, "_query_clingen_by_hgvs_batch", fake_clingen_batch)
    monkeypatch.setattr(
        mv,
        "_extract_hgvs_from_clingen",
        lambda data, transcript_nm: (data["hgvs"], f"{transcript_nm}:c.1A>G", "NP_000001.1:p.Ala1Val"),
    )
    monkeypatch.setattr(mv, "_extract_clingen_allele_id", lambda data: data.get("id"))
    monkeypatch.setattr(mv, "_clingen_allele_type", lambda data: "SNV")

    mv.map_variants(str(input_path), str(output_path))

    out_rows = _read_tsv(output_path)
    assert [r["variant_urn"] for r in out_rows] == ["v0", "v1", "v2", "v3"]


def test_map_variants_retries_on_137_with_chunking(tmp_path, monkeypatch):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"

    rows = [
        {"variant_urn": "v0", "raw_hgvs_nt": "", "raw_hgvs_pro": "p.Ala1Val", "target_sequence": "SEQ_A"},
        {"variant_urn": "v1", "raw_hgvs_nt": "", "raw_hgvs_pro": "p.Ala2Val", "target_sequence": "SEQ_A"},
        {"variant_urn": "v2", "raw_hgvs_nt": "", "raw_hgvs_pro": "p.Ala3Val", "target_sequence": "SEQ_A"},
    ]
    _write_tsv(input_path, rows)

    calls = []

    async def fake_pipeline(group_name, target_seq, row_entries, dcd):
        calls.append((group_name, len(row_entries)))
        if "#retry" not in group_name:
            raise RuntimeError("BLAT process returned error code 137")
        per_row = [(orig_idx, f"NC_000001.11:g.{orig_idx + 100}A>G", None) for orig_idx, *_ in row_entries]
        return per_row, "NM_000001.1"

    async def fake_clingen_batch(hgvs_strings, max_concurrency=5):
        return {h: {"hgvs": h, "id": "CA123"} for h in hgvs_strings}

    monkeypatch.setattr(mv, "_try_import_dcd_mapping", lambda: object())
    monkeypatch.setattr(mv, "_run_dcd_mapping_pipeline", fake_pipeline)
    monkeypatch.setattr(mv, "_query_clingen_by_hgvs_batch", fake_clingen_batch)
    monkeypatch.setattr(mv, "_extract_hgvs_from_clingen", lambda data, transcript_nm: (data["hgvs"], "NM_1:c.1A>G", "NP_1:p.Ala1Val"))
    monkeypatch.setattr(mv, "_extract_clingen_allele_id", lambda data: data.get("id"))
    monkeypatch.setattr(mv, "_clingen_allele_type", lambda data: "SNV")

    mv.map_variants(
        str(input_path),
        str(output_path),
        dcd_chunk_size_on_137=2,
        dcd_max_retry_attempts=3,
        preserve_order="index",
    )

    # First call should be full group, then retry chunks of size 2 and 1.
    assert calls[0][1] == 3
    chunk_sizes = [size for name, size in calls[1:] if "#retry1_chunk" in name]
    assert chunk_sizes == [2, 1]

    out_rows = _read_tsv(output_path)
    assert all((r.get("mapping_error") or "") == "" for r in out_rows)


def test_map_variants_no_retry_when_137_retry_disabled(tmp_path, monkeypatch):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"

    rows = [
        {"variant_urn": "v0", "raw_hgvs_nt": "", "raw_hgvs_pro": "p.Ala1Val", "target_sequence": "SEQ_A"},
    ]
    _write_tsv(input_path, rows)

    async def always_137(group_name, target_seq, row_entries, dcd):
        raise RuntimeError("BLAT process returned error code 137")

    monkeypatch.setattr(mv, "_try_import_dcd_mapping", lambda: object())
    monkeypatch.setattr(mv, "_run_dcd_mapping_pipeline", always_137)

    with pytest.raises(RuntimeError, match="error code 137"):
        mv.map_variants(
            str(input_path),
            str(output_path),
            dcd_chunk_on_137=False,
        )


@pytest.mark.parametrize(
    "raw_nt,raw_pro,expected",
    [
        ("NM_000001.1:c.123A>G", "", 1),
        ("ENST00000316054.9:c.1142G>A", "", 1),
        ("c.123A>G", "", 2),
        ("", "p.Ala1Val", 3),
        ("", "", None),
        ("_wt", "", None),
    ],
)
def test_detect_case_variants(raw_nt, raw_pro, expected):
    assert mv._detect_case(raw_nt, raw_pro) == expected


def test_map_variants_routes_class1_class2_class3(tmp_path, monkeypatch):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"

    rows = [
        {
            "variant_urn": "v_class1",
            "raw_hgvs_nt": "NM_000001.1:c.10A>G",
            "raw_hgvs_pro": "",
            "target_sequence": "SEQ_SHARED",
        },
        {
            "variant_urn": "v_class2",
            "raw_hgvs_nt": "c.20A>G",
            "raw_hgvs_pro": "",
            "target_sequence": "SEQ_SHARED",
        },
        {
            "variant_urn": "v_class3",
            "raw_hgvs_nt": "",
            "raw_hgvs_pro": "p.Ala3Val",
            "target_sequence": "SEQ_SHARED",
        },
    ]
    _write_tsv(input_path, rows)

    called_case1 = []
    pipeline_calls = []

    def fake_case1(raw_hgvs_nt, dcd):
        called_case1.append(raw_hgvs_nt)
        return (
            "NM_000001.1:c.10A>G",
            "NC_000001.11:g.10A>G",
            "NP_000001.1:p.Lys4Arg",
            None,
            "CA001",
        )

    async def fake_pipeline(group_name, target_seq, row_entries, dcd):
        pipeline_calls.append((group_name, target_seq, row_entries))
        per_row = []
        for orig_idx, hgvs_nt, hgvs_pro, case in row_entries:
            per_row.append((orig_idx, f"NC_000001.11:g.{orig_idx + 100}A>G", None))
        return per_row, "NM_000001.1"

    async def fake_clingen_batch(hgvs_strings, max_concurrency=5):
        return {h: {"hgvs": h, "id": "CA777"} for h in hgvs_strings}

    monkeypatch.setattr(mv, "_try_import_dcd_mapping", lambda: object())
    monkeypatch.setattr(mv, "_process_case1", fake_case1)
    monkeypatch.setattr(mv, "_run_dcd_mapping_pipeline", fake_pipeline)
    monkeypatch.setattr(mv, "_query_clingen_by_hgvs_batch", fake_clingen_batch)
    monkeypatch.setattr(mv, "_extract_hgvs_from_clingen", lambda data, transcript_nm: (data["hgvs"], "NM_000001.1:c.99A>G", "NP_000001.1:p.Arg33Gly"))
    monkeypatch.setattr(mv, "_extract_clingen_allele_id", lambda data: data.get("id"))
    monkeypatch.setattr(mv, "_clingen_allele_type", lambda data: "SNV")

    mv.map_variants(str(input_path), str(output_path), max_clingen_concurrency=1)

    assert called_case1 == ["NM_000001.1:c.10A>G"]
    assert len(pipeline_calls) == 1
    _, _, row_entries = pipeline_calls[0]
    assert sorted(case for _, _, _, case in row_entries) == [2, 3]

    out_rows = _read_tsv(output_path)
    assert [r["variant_urn"] for r in out_rows] == ["v_class1", "v_class2", "v_class3"]
    assert out_rows[0]["clingen_allele_id"] == "CA001"
    assert out_rows[1]["mapped_hgvs_g"].startswith("NC_000001.11:g.")
    assert out_rows[2]["mapped_hgvs_g"].startswith("NC_000001.11:g.")


def test_map_variants_merge_existing_reuses_rows(tmp_path, monkeypatch):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"
    existing_path = tmp_path / "existing.tsv"

    input_rows = [
        {
            "variant_urn": "v0",
            "raw_hgvs_nt": "",
            "raw_hgvs_pro": "p.Ala1Val",
            "target_sequence": "SEQ_A",
        }
    ]
    _write_tsv(input_path, input_rows)

    with open(existing_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "variant_urn",
                "raw_hgvs_nt",
                "raw_hgvs_pro",
                "target_sequence",
                "mapped_hgvs_g",
                "mapped_hgvs_c",
                "mapped_hgvs_p",
                "mapping_error",
                "clingen_allele_id",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {
                "variant_urn": "v0",
                "raw_hgvs_nt": "",
                "raw_hgvs_pro": "p.Ala1Val",
                "target_sequence": "WHATEVER",
                "mapped_hgvs_g": "NC_000001.11:g.123A>G",
                "mapped_hgvs_c": "NM_000001.1:c.123A>G",
                "mapped_hgvs_p": "NP_000001.1:p.Ala41Val",
                "mapping_error": "",
                "clingen_allele_id": "CA999",
            }
        )

    async def should_not_run(*args, **kwargs):
        raise AssertionError("Pipeline should not run for merged rows")

    monkeypatch.setattr(mv, "_run_dcd_mapping_pipeline", should_not_run)
    monkeypatch.setattr(mv, "_try_import_dcd_mapping", lambda: object())

    mv.map_variants(
        str(input_path),
        str(output_path),
        merge_existing_files=(str(existing_path),),
    )

    out_rows = _read_tsv(output_path)
    assert len(out_rows) == 1
    assert out_rows[0]["clingen_allele_id"] == "CA999"
    assert out_rows[0]["mapped_hgvs_g"] == "NC_000001.11:g.123A>G"


def test_process_case1_batch_deduplicates_assays_and_preserves_order(monkeypatch):
    queried = []

    dcd = {
        "fetch_clingen_genomic_hgvs": lambda raw: "NC_000001.11:g.111A>G" if "c.10" in raw else "NC_000001.11:g.222A>T"
    }

    async def fake_query_batch(assays, max_concurrency=5):
        queried.extend(assays)
        return {assay: {"assay": assay, "id": f"CA-{assay.split(':')[-1]}"} for assay in assays}

    monkeypatch.setattr(mv, "_query_clingen_by_hgvs_batch", fake_query_batch)
    monkeypatch.setattr(
        mv,
        "_extract_hgvs_from_clingen",
        lambda data, original_accession: (data["assay"], f"{original_accession}:c.1A>G", "NP_000001.1:p.Lys1Arg"),
    )
    monkeypatch.setattr(mv, "_extract_clingen_allele_id", lambda data: data["id"])

    raws = [
        "NM_000001.1:c.10A>G",  # maps to assay 111
        "NM_000002.1:c.20A>T",  # maps to assay 222
        "NM_000003.1:c.10A>G",  # maps again to assay 111 (dedupe expected)
    ]
    results = mv._process_case1_batch(raws, dcd=dcd, max_concurrency=3)

    assert queried == ["NC_000001.11:g.111A>G", "NC_000001.11:g.222A>T"]
    assert len(results) == 3
    assert [r[1] for r in results] == [
        "NC_000001.11:g.111A>G",
        "NC_000001.11:g.222A>T",
        "NC_000001.11:g.111A>G",
    ]
    assert all(r[3] is None for r in results)


def test_process_case1_batch_reports_invalid_and_missing_clingen(monkeypatch):
    async def fake_query_batch(assays, max_concurrency=5):
        # Return no data for the queried assay to exercise missing ClinGen handling.
        return {assays[0]: None}

    monkeypatch.setattr(mv, "_query_clingen_by_hgvs_batch", fake_query_batch)

    raws = [
        "bad-format",
        "NM_000001.1:c.10A>G",
    ]
    results = mv._process_case1_batch(raws, dcd=None, max_concurrency=2)

    assert "Expected 'accession:variant' format" in (results[0][3] or "")
    assert "ClinGen returned no data" in (results[1][3] or "")
