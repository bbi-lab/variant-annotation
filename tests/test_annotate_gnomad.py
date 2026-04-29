import csv
import time

import pytest
import src.annotate_gnomad as mod

from src.annotate_gnomad import GnomadRecord, annotate_row


def test_annotate_row_emits_pipe_aligned_values_for_all_candidates():
    records = {
        "CA2": GnomadRecord(
            caid="CA2",
            allele_count=10,
            allele_number=100,
            allele_frequency=0.1,
            minor_allele_frequency=0.1,
            faf95_max=0.09,
            faf95_max_ancestry="nfe",
        ),
        "CA3": GnomadRecord(
            caid="CA3",
            allele_count=5,
            allele_number=100,
            allele_frequency=0.05,
            minor_allele_frequency=0.05,
            faf95_max=None,
            faf95_max_ancestry="",
        ),
    }

    row = {"dna_clingen_allele_id": "CA1|CA2|CA3"}
    out = annotate_row(row, records, "gnomad.v4.1", "dna_clingen_allele_id")

    assert out["gnomad.v4.1.minor_allele_frequency"] == "|0.1|0.05"
    assert out["gnomad.v4.1.allele_frequency"] == "|0.1|0.05"
    assert out["gnomad.v4.1.allele_count"] == "|10|5"
    assert out["gnomad.v4.1.allele_number"] == "|100|100"
    assert out["gnomad.v4.1.faf95_max"] == "|0.09|"
    assert out["gnomad.v4.1.faf95_max_ancestry"] == "|nfe|"


def test_annotate_row_handles_no_match():
    records = {}
    row = {"dna_clingen_allele_id": "CA1|CA2"}
    out = annotate_row(row, records, "gnomad.v4.1", "dna_clingen_allele_id")

    assert out["gnomad.v4.1.minor_allele_frequency"] == "|"
    assert out["gnomad.v4.1.allele_frequency"] == "|"


def test_annotate_row_custom_dna_column():
    records = {
        "CA7": GnomadRecord(
            caid="CA7",
            allele_count=2,
            allele_number=20,
            allele_frequency=0.1,
            minor_allele_frequency=0.1,
            faf95_max=None,
            faf95_max_ancestry="",
        )
    }
    row = {"my_dna_ids": "CA7"}
    out = annotate_row(row, records, "gnomad.v4.1", "my_dna_ids")

    assert out["gnomad.v4.1.minor_allele_frequency"] == "0.1"


def test_annotate_row_preserves_empty_candidate_slots():
    records = {
        "CA7": GnomadRecord(
            caid="CA7",
            allele_count=2,
            allele_number=20,
            allele_frequency=0.1,
            minor_allele_frequency=0.1,
            faf95_max=None,
            faf95_max_ancestry="",
        )
    }

    row = {"dna_clingen_allele_id": "CA7||CA7"}
    out = annotate_row(row, records, "gnomad.v4.1", "dna_clingen_allele_id")

    assert out["gnomad.v4.1.minor_allele_frequency"] == "0.1||0.1"


def test_main_applies_skip_and_limit(tmp_path, monkeypatch):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"

    with in_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["variant_urn", "dna_clingen_allele_id"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow({"variant_urn": "v1", "dna_clingen_allele_id": "CA1"})
        writer.writerow({"variant_urn": "v2", "dna_clingen_allele_id": "CA2"})
        writer.writerow({"variant_urn": "v3", "dna_clingen_allele_id": "CA3"})
        writer.writerow({"variant_urn": "v4", "dna_clingen_allele_id": "CA4"})

    monkeypatch.setattr(mod, "ensure_local_gnomad_ht", lambda *args, **kwargs: tmp_path / "dummy.ht")
    monkeypatch.setattr(
        mod,
        "load_gnomad_records_for_caids",
        lambda local_ht, caids, cache_dir, **kwargs: {
            caid: GnomadRecord(
                caid=caid,
                allele_count=1,
                allele_number=10,
                allele_frequency=0.1,
                minor_allele_frequency=0.1,
                faf95_max=None,
                faf95_max_ancestry="",
            )
            for caid in caids
        },
    )

    mod.main([
        str(in_path),
        str(out_path),
        "--gnomad-ht-uri",
        "dummy-uri",
        "--skip",
        "1",
        "--limit",
        "2",
    ])

    with out_path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    assert [r["variant_urn"] for r in rows] == ["v2", "v3"]
    assert [r["gnomad.v4_1.minor_allele_frequency"] for r in rows] == ["0.1", "0.1"]


def test_raise_actionable_hail_error_for_missing_gs_filesystem():
    exc = RuntimeError('UnsupportedFileSystemException: No FileSystem for scheme "gs"')
    with pytest.raises(RuntimeError, match="Google Cloud Storage connector"):
        mod._raise_actionable_hail_error("gs://bucket/table.ht", exc)


def test_load_gnomad_records_resolves_caids_via_clingen(monkeypatch, tmp_path):
    """load_gnomad_records_for_caids resolves each CAID to a gnomad_key via ClinGen,
    then looks up matching rows in the local Hail table."""
    import src.lib.clingen as clingen_mod

    coord_responses = {
        "CA1": ("chr7", 117548628, "G", "A"),
        "CA2": None,  # no GRCh38 coords available
    }

    def fake_resolve(caid, cache, **kwargs):
        result = coord_responses.get(caid)
        cache[caid] = result
        return result

    monkeypatch.setattr(clingen_mod, "resolve_grch38_coordinates", fake_resolve)

    # Stub out Hail: replace _import_hail with a minimal fake
    class _FakeRow:
        def __init__(self, gnomad_key, ac, an, faf, faf_anc):
            self.gnomad_key = gnomad_key
            self.allele_count = ac
            self.allele_number = an
            self.faf95_max = faf
            self.faf95_max_ancestry = faf_anc

    _rows = [_FakeRow("chr7:117548628:G:A", 20, 200, 0.08, "nfe")]

    class _FakeExpr:
        def contains(self, _other):
            return True

    class _FakeHt:
        def __init__(self):
            self.gnomad_key = None  # attribute access for expression building

        def filter(self, _expr):
            return self

        def collect(self):
            return _rows

    class _FakeHl:
        def __init__(self):
            self._ht = _FakeHt()

        def init(self, **kwargs):
            pass

        def stop(self):
            pass

        def read_table(self, path):
            return self._ht

        def literal(self, val):
            return _FakeExpr()

    fake_hl = _FakeHl()
    monkeypatch.setattr(mod, "_import_hail", lambda: fake_hl)
    monkeypatch.setattr(mod, "_hail_init_kwargs", lambda *a, **kw: {})

    result = mod.load_gnomad_records_for_caids(tmp_path / "dummy.ht", {"CA1", "CA2"}, tmp_path)

    assert "CA1" in result
    rec = result["CA1"]
    assert rec.allele_count == 20
    assert rec.allele_number == 200
    assert abs(rec.allele_frequency - 0.1) < 1e-9
    assert rec.faf95_max == 0.08
    assert rec.faf95_max_ancestry == "nfe"
    assert "CA2" not in result  # no GRCh38 coords → not resolved


def test_build_caid_to_gnomad_key_uses_input_coordinate_columns():
    rows = [
        {
            "dna_clingen_allele_id": "CA1|CA2",
            "mapped_hgvs_g_chromosome": "7|13",
            "mapped_hgvs_g_stop": "117548628|32316461",
            "mapped_hgvs_g_ref": "G|C",
            "mapped_hgvs_g_alt": "A|T",
        }
    ]
    out = mod._build_caid_to_gnomad_key(
        rows,
        dna_col="dna_clingen_allele_id",
        chrom_col="mapped_hgvs_g_chromosome",
        pos_col="mapped_hgvs_g_stop",
        ref_col="mapped_hgvs_g_ref",
        alt_col="mapped_hgvs_g_alt",
    )
    assert out["CA1"] == "chr7:117548628:G:A"
    assert out["CA2"] == "chr13:32316461:C:T"


def test_main_athena_mode_requires_output_location(tmp_path):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    with in_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["variant_urn", "dna_clingen_allele_id"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow({"variant_urn": "v1", "dna_clingen_allele_id": "CA1"})

    with pytest.raises(SystemExit):
        mod.main([
            str(in_path),
            str(out_path),
            "--execution-mode",
            "athena",
            "--gnomad-ht-uri",
            "dummy-uri",
        ])


def test_main_athena_mode_uses_athena_loader(tmp_path, monkeypatch):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"

    with in_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["variant_urn", "dna_clingen_allele_id"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow({"variant_urn": "v1", "dna_clingen_allele_id": "CA1"})

    monkeypatch.setattr(mod, "ensure_local_gnomad_ht", lambda *args, **kwargs: (_ for _ in ()).throw(AssertionError("hail path should not be used")))

    def _fake_athena_loader(caids, **kwargs):
        assert caids == {"CA1"}
        return {
            "CA1": GnomadRecord(
                caid="CA1",
                allele_count=2,
                allele_number=10,
                allele_frequency=0.2,
                minor_allele_frequency=0.2,
                faf95_max=0.1,
                faf95_max_ancestry="nfe",
            )
        }

    monkeypatch.setattr(mod, "load_gnomad_records_for_caids_athena", _fake_athena_loader)

    mod.main([
        str(in_path),
        str(out_path),
        "--execution-mode",
        "athena",
        "--athena-output-location",
        "s3://dummy-athena-output/",
        "--gnomad-ht-uri",
        "dummy-uri",
    ])

    with out_path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    assert len(rows) == 1
    assert rows[0]["gnomad.v4_1.minor_allele_frequency"] == "0.2"


def test_main_athena_mode_streams_batches_in_input_order(tmp_path, monkeypatch):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"

    with in_path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["variant_urn", "dna_clingen_allele_id"],
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerow({"variant_urn": "v1", "dna_clingen_allele_id": "CA1"})
        writer.writerow({"variant_urn": "v2", "dna_clingen_allele_id": "CA2"})
        writer.writerow({"variant_urn": "v3", "dna_clingen_allele_id": "CA1"})

    calls: list[set[str]] = []

    def _fake_athena_loader(caids, **kwargs):
        calls.append(set(caids))
        out = {}
        for caid in caids:
            if caid == "CA1":
                out[caid] = GnomadRecord(
                    caid="CA1",
                    allele_count=2,
                    allele_number=10,
                    allele_frequency=0.2,
                    minor_allele_frequency=0.2,
                    faf95_max=0.1,
                    faf95_max_ancestry="nfe",
                )
            elif caid == "CA2":
                out[caid] = GnomadRecord(
                    caid="CA2",
                    allele_count=1,
                    allele_number=10,
                    allele_frequency=0.1,
                    minor_allele_frequency=0.1,
                    faf95_max=0.05,
                    faf95_max_ancestry="afr",
                )
        return out

    monkeypatch.setattr(mod, "load_gnomad_records_for_caids_athena", _fake_athena_loader)

    mod.main([
        str(in_path),
        str(out_path),
        "--execution-mode",
        "athena",
        "--athena-output-location",
        "s3://dummy-athena-output/",
        "--athena-row-batch-size",
        "1",
        "--gnomad-ht-uri",
        "dummy-uri",
    ])

    with out_path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    assert [r["variant_urn"] for r in rows] == ["v1", "v2", "v3"]
    assert [r["gnomad.v4_1.minor_allele_frequency"] for r in rows] == ["0.2", "0.1", "0.2"]
    # CA1 should only be looked up once across batches thanks to cache reuse.
    assert calls == [{"CA1"}, {"CA2"}]


def test_cache_progress_message_reports_file_growth(tmp_path):
    cache_dir = tmp_path / "gnomad_cache"
    ht_path = cache_dir / "gnomad_v4_1_indexed.ht"
    ht_path.mkdir(parents=True)
    (ht_path / "part-00000").write_text("abc", encoding="utf-8")

    message = mod._cache_progress_message("writing local cache table", time.monotonic() - 30, cache_dir, ht_path)

    assert "stage=writing local cache table" in message
    assert "cache_dir=" in message
    assert "local_ht=" in message
    assert "success_marker=no" in message
