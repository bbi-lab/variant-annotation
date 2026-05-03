import csv
import asyncio

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


def test_normalize_case2_haplotype_to_delins_success(monkeypatch):
    async def fake_run_pipeline(group_name, target_sequence, row_entries, dcd, allow_row_fallback=True, precomputed_align_result=None, precomputed_transcript=None):
        # Component mapping call from haplotype normalization helper.
        assert group_name.endswith("#haplotype-components")
        assert target_sequence == "ATGCC"
        return [(0, "NC_1", ""), (1, "NC_2", "")], "NM_000001.1"

    def fake_query_clingen(hgvs):
        if hgvs == "NC_1":
            return {"c": "NM_000001.1:c.1A>G"}
        if hgvs == "NC_2":
            return {"c": "NM_000001.1:c.3G>T"}
        return None

    monkeypatch.setattr(mv, "_run_dcd_mapping_pipeline", fake_run_pipeline)
    monkeypatch.setattr(mv, "_query_clingen_by_hgvs", fake_query_clingen)
    monkeypatch.setattr(mv, "_extract_hgvs_from_clingen", lambda data, tx: (None, data.get("c"), None))

    normalized, error = asyncio.run(
        mv._normalize_case2_haplotype_to_delins(
            raw_hgvs_nt="c.[1A>G;3G>T]",
            target_sequence="ATGCC",
            row_label="row-1",
            dcd={},
        )
    )

    assert error is None
    assert normalized == "c.1_3delinsGTT"


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


def test_map_variants_groups_mode_uses_contiguous_blocks(tmp_path, monkeypatch):
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"

    rows = [
        {"variant_urn": "v0", "raw_hgvs_nt": "", "raw_hgvs_pro": "p.Ala1Val", "target_sequence": "SEQ_A"},
        {"variant_urn": "v1", "raw_hgvs_nt": "", "raw_hgvs_pro": "p.Ala2Val", "target_sequence": "SEQ_B"},
        {"variant_urn": "v2", "raw_hgvs_nt": "", "raw_hgvs_pro": "p.Ala3Val", "target_sequence": "SEQ_A"},
    ]
    _write_tsv(input_path, rows)

    calls = []

    async def fake_pipeline(group_name, target_seq, row_entries, dcd):
        calls.append((group_name, [orig_idx for orig_idx, *_ in row_entries]))
        per_row = []
        for orig_idx, *_ in row_entries:
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

    mv.map_variants(str(input_path), str(output_path), preserve_order="groups")

    assert calls == [("SEQ_A", [0]), ("SEQ_B", [1]), ("SEQ_A", [2])]
    out_rows = _read_tsv(output_path)
    assert [r["variant_urn"] for r in out_rows] == ["v0", "v1", "v2"]


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


# ---------------------------------------------------------------------------
# Targets-file tests
# ---------------------------------------------------------------------------


def _write_targets_tsv(path, rows, fieldnames=None):
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    with open(path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_load_targets_file_basic(tmp_path):
    targets_path = tmp_path / "targets.tsv"
    _write_targets_tsv(
        targets_path,
        [
            {"target_name": "geneA", "target_sequence": "ACGT", "offset": "0"},
            {"target_name": "geneB", "target_sequence": "TTTT", "offset": "5"},
        ],
    )
    result = mv._load_targets_file(str(targets_path), "target_name")
    assert set(result) == {"geneA", "geneB"}
    assert result["geneA"]["target_sequence"] == "ACGT"
    assert result["geneB"]["offset"] == "5"


def test_load_targets_file_missing_name_col(tmp_path):
    targets_path = tmp_path / "targets.tsv"
    _write_targets_tsv(targets_path, [{"seq": "ACGT"}])
    with pytest.raises(ValueError, match="target_name"):
        mv._load_targets_file(str(targets_path), "target_name")


def test_map_variants_with_targets_file_populates_sequence(tmp_path, monkeypatch):
    """target_sequence is filled from the targets file when absent from the input."""
    targets_path = tmp_path / "targets.tsv"
    _write_targets_tsv(
        targets_path,
        [{"target_name": "geneA", "target_sequence": "SEQ_A"}],
    )

    input_path = tmp_path / "in.tsv"
    with open(input_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["variant_urn", "target_name", "raw_hgvs_nt", "raw_hgvs_pro"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow({"variant_urn": "v0", "target_name": "geneA", "raw_hgvs_nt": "", "raw_hgvs_pro": "p.Ala1Val"})

    output_path = tmp_path / "out.tsv"

    async def fake_pipeline(group_name, target_seq, row_entries, dcd):
        assert target_seq == "SEQ_A", "target_sequence was not merged from targets file"
        return [(orig_idx, f"NC_000001.11:g.{orig_idx}A>G", None) for orig_idx, *_ in row_entries], "NM_000001.1"

    async def fake_clingen_batch(hgvs_strings, max_concurrency=5):
        return {h: {"hgvs": h, "id": "CA1"} for h in hgvs_strings}

    monkeypatch.setattr(mv, "_try_import_dcd_mapping", lambda: object())
    monkeypatch.setattr(mv, "_run_dcd_mapping_pipeline", fake_pipeline)
    monkeypatch.setattr(mv, "_query_clingen_by_hgvs_batch", fake_clingen_batch)
    monkeypatch.setattr(mv, "_extract_hgvs_from_clingen", lambda data, tx: (data["hgvs"], None, None))
    monkeypatch.setattr(mv, "_extract_clingen_allele_id", lambda data: data.get("id"))
    monkeypatch.setattr(mv, "_clingen_allele_type", lambda data: "CA")

    mv.map_variants(
        str(input_path),
        str(output_path),
        targets_file=str(targets_path),
        target_name_col="target_name",
    )

    out_rows = _read_tsv(output_path)
    assert len(out_rows) == 1
    assert out_rows[0]["variant_urn"] == "v0"
    # target_sequence from targets file should appear in output
    assert out_rows[0]["target_sequence"] == "SEQ_A"


def test_map_variants_with_targets_file_extra_cols_in_output(tmp_path, monkeypatch):
    """Extra columns from the targets file appear in the output."""
    targets_path = tmp_path / "targets.tsv"
    _write_targets_tsv(
        targets_path,
        [{"target_name": "geneA", "target_sequence": "SEQ_A", "uniprot_id": "P12345"}],
    )

    input_path = tmp_path / "in.tsv"
    with open(input_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["variant_urn", "target_name", "raw_hgvs_nt", "raw_hgvs_pro"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow({"variant_urn": "v0", "target_name": "geneA", "raw_hgvs_nt": "", "raw_hgvs_pro": "p.Ala1Val"})

    output_path = tmp_path / "out.tsv"

    async def fake_pipeline(group_name, target_seq, row_entries, dcd):
        return [(orig_idx, f"NC_000001.11:g.{orig_idx}A>G", None) for orig_idx, *_ in row_entries], "NM_000001.1"

    async def fake_clingen_batch(hgvs_strings, max_concurrency=5):
        return {h: {"hgvs": h, "id": "CA1"} for h in hgvs_strings}

    monkeypatch.setattr(mv, "_try_import_dcd_mapping", lambda: object())
    monkeypatch.setattr(mv, "_run_dcd_mapping_pipeline", fake_pipeline)
    monkeypatch.setattr(mv, "_query_clingen_by_hgvs_batch", fake_clingen_batch)
    monkeypatch.setattr(mv, "_extract_hgvs_from_clingen", lambda data, tx: (data["hgvs"], None, None))
    monkeypatch.setattr(mv, "_extract_clingen_allele_id", lambda data: data.get("id"))
    monkeypatch.setattr(mv, "_clingen_allele_type", lambda data: "CA")

    mv.map_variants(
        str(input_path),
        str(output_path),
        targets_file=str(targets_path),
        target_name_col="target_name",
    )

    out_rows = _read_tsv(output_path)
    assert out_rows[0]["uniprot_id"] == "P12345"


def test_map_variants_targets_file_unknown_name_logs_warning(tmp_path, monkeypatch, caplog):
    """A row with an unrecognised target_name logs a warning and gets a blank sequence."""
    import logging

    targets_path = tmp_path / "targets.tsv"
    _write_targets_tsv(targets_path, [{"target_name": "geneA", "target_sequence": "SEQ_A"}])

    input_path = tmp_path / "in.tsv"
    with open(input_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["variant_urn", "target_name", "raw_hgvs_nt", "raw_hgvs_pro"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {"variant_urn": "v0", "target_name": "UNKNOWN", "raw_hgvs_nt": "", "raw_hgvs_pro": "p.Ala1Val"}
        )

    output_path = tmp_path / "out.tsv"

    with caplog.at_level(logging.WARNING, logger="src.map_variants"):
        mv.map_variants(
            str(input_path),
            str(output_path),
            targets_file=str(targets_path),
            target_name_col="target_name",
        )

    assert any("UNKNOWN" in record.message for record in caplog.records)


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


# ---------------------------------------------------------------------------
# normalize_protein_hgvs tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "input_hgvs,expected",
    [
        # Substitution
        ("p.A300T", "p.Ala300Thr"),
        ("p.K57N", "p.Lys57Asn"),
        ("p.M1V", "p.Met1Val"),
        # Synonymous (= → ref AA repeated)
        ("p.A300=", "p.Ala300Ala"),
        ("p.K57=", "p.Lys57Lys"),
        # Stop gain (* → Ter)
        ("p.A300*", "p.Ala300Ter"),
        ("p.Q192*", "p.Gln192Ter"),
        # Deletion via dash
        ("p.A300-", "p.Ala300del"),
        # Deletion via keyword – ref AA still converted
        ("p.A300del", "p.Ala300del"),
        # Frameshift – suffix kept verbatim
        ("p.A300fs", "p.Ala300fs"),
        ("p.A300fs*7", "p.Ala300fs*7"),
        # Insertion keyword
        ("p.A300ins", "p.Ala300ins"),
        # Duplication keyword
        ("p.A300dup", "p.Ala300dup"),
        # Already 3-letter – returned unchanged
        ("p.Ala300Thr", "p.Ala300Thr"),
        ("p.Lys57Asn", "p.Lys57Asn"),
        ("p.Ala300del", "p.Ala300del"),
        # Prefixed with accession
        ("NP_000001.1:p.A300T", "NP_000001.1:p.Ala300Thr"),
        # Blank / sentinel – returned as-is
        ("", ""),
        ("_wt", "_wt"),
        # Non-matching pattern – returned as-is
        ("c.300A>T", "c.300A>T"),
    ],
)
def test_normalize_protein_hgvs(input_hgvs, expected):
    assert mv.normalize_protein_hgvs(input_hgvs) == expected


def test_map_variants_normalizes_raw_hgvs_pro_in_output(tmp_path, monkeypatch):
    """1-letter p. strings are normalized to 3-letter in the output row."""
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"

    with open(input_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["variant_urn", "raw_hgvs_nt", "raw_hgvs_pro", "target_sequence"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {"variant_urn": "v0", "raw_hgvs_nt": "", "raw_hgvs_pro": "p.A300T", "target_sequence": "SEQ_A"}
        )

    async def fake_pipeline(group_name, target_seq, row_entries, dcd):
        return [(orig_idx, f"NC_000001.11:g.{orig_idx}A>G", None) for orig_idx, *_ in row_entries], "NM_1"

    async def fake_clingen_batch(hgvs_strings, max_concurrency=5):
        return {h: {"hgvs": h, "id": "CA1"} for h in hgvs_strings}

    monkeypatch.setattr(mv, "_try_import_dcd_mapping", lambda: object())
    monkeypatch.setattr(mv, "_run_dcd_mapping_pipeline", fake_pipeline)
    monkeypatch.setattr(mv, "_query_clingen_by_hgvs_batch", fake_clingen_batch)
    monkeypatch.setattr(mv, "_extract_hgvs_from_clingen", lambda data, tx: (data["hgvs"], None, "p.Ala300Thr"))
    monkeypatch.setattr(mv, "_extract_clingen_allele_id", lambda data: data.get("id"))
    monkeypatch.setattr(mv, "_clingen_allele_type", lambda data: "CA")

    mv.map_variants(str(input_path), str(output_path))

    out_rows = _read_tsv(output_path)
    assert out_rows[0]["raw_hgvs_pro"] == "p.Ala300Thr"


# ---------------------------------------------------------------------------
# _extract_clingen_allele_id tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "data,expected",
    [
        # Normal CA / PA IDs
        ({"@id": "https://reg.genome.network/allele/CA123456"}, "CA123456"),
        ({"@id": "https://reg.genome.network/allele/PA987654"}, "PA987654"),
        # Blank-node placeholders must be rejected
        ({"@id": "_:PA3262009431"}, None),
        ({"@id": "_:CA0000001"}, None),
        # Fallback via "id" field
        ({"id": "CA123456"}, "CA123456"),
        ({"id": "_:PA123"}, None),
        # Empty / missing
        ({"@id": ""}, None),
        ({}, None),
    ],
)
def test_extract_clingen_allele_id(data, expected):
    assert mv._extract_clingen_allele_id(data) == expected


def test_map_variants_clingen_no_data_preserves_protein_hgvs(tmp_path, monkeypatch):
    """When ClinGen returns no data, the dcd_mapping p. result is kept in mapped_hgvs_p."""
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"

    with open(input_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["variant_urn", "raw_hgvs_nt", "raw_hgvs_pro", "target_sequence"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {"variant_urn": "v0", "raw_hgvs_nt": "", "raw_hgvs_pro": "p.Q2Ter", "target_sequence": "SEQ_A"}
        )

    async def fake_pipeline(group_name, target_seq, row_entries, dcd):
        return [(orig_idx, "NP_005624.2:p.Gln2Ter", None) for orig_idx, *_ in row_entries], "NM_005633.4"

    async def fake_clingen_batch(hgvs_strings, max_concurrency=5):
        # Simulate ClinGen returning no data
        return {h: None for h in hgvs_strings}

    monkeypatch.setattr(mv, "_try_import_dcd_mapping", lambda: object())
    monkeypatch.setattr(mv, "_run_dcd_mapping_pipeline", fake_pipeline)
    monkeypatch.setattr(mv, "_query_clingen_by_hgvs_batch", fake_clingen_batch)

    mv.map_variants(str(input_path), str(output_path))

    out_rows = _read_tsv(output_path)
    assert out_rows[0]["mapped_hgvs_p"] == "NP_005624.2:p.Gln2Ter"
    assert "ClinGen returned no data" in out_rows[0]["mapping_error"]


def test_map_variants_clingen_no_data_preserves_genomic_hgvs(tmp_path, monkeypatch):
    """When ClinGen returns no data, the dcd_mapping genomic result is kept in mapped_hgvs_g."""
    input_path = tmp_path / "in.tsv"
    output_path = tmp_path / "out.tsv"

    with open(input_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["variant_urn", "raw_hgvs_nt", "raw_hgvs_pro", "target_sequence"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerow(
            {"variant_urn": "v0", "raw_hgvs_nt": "", "raw_hgvs_pro": "p.Q2Ter", "target_sequence": "SEQ_A"}
        )

    genomic = "NC_000002.12:g.12345A>T"

    async def fake_pipeline(group_name, target_seq, row_entries, dcd):
        return [(orig_idx, genomic, None) for orig_idx, *_ in row_entries], "NM_005633.4"

    async def fake_clingen_batch(hgvs_strings, max_concurrency=5):
        return {h: None for h in hgvs_strings}

    monkeypatch.setattr(mv, "_try_import_dcd_mapping", lambda: object())
    monkeypatch.setattr(mv, "_run_dcd_mapping_pipeline", fake_pipeline)
    monkeypatch.setattr(mv, "_query_clingen_by_hgvs_batch", fake_clingen_batch)

    mv.map_variants(str(input_path), str(output_path))

    out_rows = _read_tsv(output_path)
    assert out_rows[0]["mapped_hgvs_g"] == genomic
    assert "ClinGen returned no data" in out_rows[0]["mapping_error"]


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
