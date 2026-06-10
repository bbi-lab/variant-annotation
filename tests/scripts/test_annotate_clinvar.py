"""Unit tests for src/annotate_clinvar.py."""

from __future__ import annotations

import csv
import gzip
import io
from pathlib import Path
from unittest.mock import patch

import pytest
import annotate_clinvar as mod

from annotate_clinvar import (
    annotate_row,
    load_clinvar_tsv,
    stars_for_review_status,
)



# ---------------------------------------------------------------------------
# stars_for_review_status
# ---------------------------------------------------------------------------


@pytest.mark.unit
class TestStarsForReviewStatus:
    def test_practice_guideline_4(self):
        assert stars_for_review_status("practice guideline") == "4"

    def test_expert_panel_3(self):
        assert stars_for_review_status("reviewed by expert panel") == "3"

    def test_multiple_submitters_no_conflicts_2(self):
        assert stars_for_review_status("criteria provided, multiple submitters, no conflicts") == "2"

    def test_single_submitter_1(self):
        assert stars_for_review_status("criteria provided, single submitter") == "1"

    def test_conflicting_classifications_1(self):
        assert stars_for_review_status("criteria provided, conflicting classifications") == "1"

    def test_no_assertion_criteria_0(self):
        assert stars_for_review_status("no assertion criteria provided") == "0"

    def test_no_assertion_0(self):
        assert stars_for_review_status("no assertion provided") == "0"

    def test_empty_string_returns_empty(self):
        assert stars_for_review_status("") == ""

    def test_unknown_value_returns_empty(self):
        assert stars_for_review_status("some unknown review status") == ""

    def test_case_insensitive(self):
        assert stars_for_review_status("Practice Guideline") == "4"

    def test_leading_trailing_whitespace(self):
        assert stars_for_review_status("  reviewed by expert panel  ") == "3"


# ---------------------------------------------------------------------------
# load_clinvar_tsv
# ---------------------------------------------------------------------------


def _make_gz_tsv(rows: list[dict], fieldnames: list[str]) -> bytes:
    """Create an in-memory gzipped TSV from rows."""
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    for row in rows:
        writer.writerow(row)
    data = buf.getvalue().encode("utf-8")
    gz_buf = io.BytesIO()
    with gzip.GzipFile(fileobj=gz_buf, mode="wb") as gz:
        gz.write(data)
    return gz_buf.getvalue()


FIELDNAMES = [
    "#AlleleID",
    "ClinicalSignificance",
    "ReviewStatus",
    "LastEvaluated",
    "Assembly",
    "Origin",
]


@pytest.mark.unit
class TestLoadClinvarTsv:
    def test_loads_basic_record(self, tmp_path):
        gz_bytes = _make_gz_tsv(
            [
                {
                    "#AlleleID": "12345",
                    "ClinicalSignificance": "Pathogenic",
                    "ReviewStatus": "reviewed by expert panel",
                    "LastEvaluated": "2023-01-15",
                    "Assembly": "GRCh38",
                    "Origin": "germline",
                }
            ],
            FIELDNAMES,
        )
        path = tmp_path / "clinvar.txt.gz"
        path.write_bytes(gz_bytes)
        data = load_clinvar_tsv(path)
        assert "12345" in data
        assert data["12345"]["ClinicalSignificance"] == "Pathogenic"
        assert data["12345"]["ReviewStatus"] == "reviewed by expert panel"
        assert data["12345"]["LastEvaluated"] == "2023-01-15"

    def test_skips_non_grch38_rows(self, tmp_path):
        gz_bytes = _make_gz_tsv(
            [
                {
                    "#AlleleID": "111",
                    "ClinicalSignificance": "Benign",
                    "ReviewStatus": "criteria provided, single submitter",
                    "LastEvaluated": "2022-06-01",
                    "Assembly": "GRCh37",
                    "Origin": "germline",
                },
                {
                    "#AlleleID": "222",
                    "ClinicalSignificance": "Pathogenic",
                    "ReviewStatus": "reviewed by expert panel",
                    "LastEvaluated": "2023-03-01",
                    "Assembly": "GRCh38",
                    "Origin": "germline",
                },
            ],
            FIELDNAMES,
        )
        path = tmp_path / "clinvar.txt.gz"
        path.write_bytes(gz_bytes)
        data = load_clinvar_tsv(path)
        assert "111" not in data
        assert "222" in data

    def test_empty_allele_id_skipped(self, tmp_path):
        gz_bytes = _make_gz_tsv(
            [
                {
                    "#AlleleID": "",
                    "ClinicalSignificance": "Pathogenic",
                    "ReviewStatus": "no assertion criteria provided",
                    "LastEvaluated": "",
                    "Assembly": "GRCh38",
                    "Origin": "germline",
                }
            ],
            FIELDNAMES,
        )
        path = tmp_path / "clinvar.txt.gz"
        path.write_bytes(gz_bytes)
        data = load_clinvar_tsv(path)
        assert data == {}

    def test_multiple_records(self, tmp_path):
        rows = [
            {
                "#AlleleID": str(i),
                "ClinicalSignificance": "Pathogenic",
                "ReviewStatus": "criteria provided, single submitter",
                "LastEvaluated": "2024-01-01",
                "Assembly": "GRCh38",
                "Origin": "germline",
            }
            for i in range(1, 6)
        ]
        gz_bytes = _make_gz_tsv(rows, FIELDNAMES)
        path = tmp_path / "clinvar.txt.gz"
        path.write_bytes(gz_bytes)
        data = load_clinvar_tsv(path)
        assert len(data) == 5

    def test_normalises_dash_placeholders_to_empty(self, tmp_path):
        gz_bytes = _make_gz_tsv(
            [
                {
                    "#AlleleID": "333",
                    "ClinicalSignificance": "-",
                    "ReviewStatus": "-",
                    "LastEvaluated": "-",
                    "Assembly": "GRCh38",
                    "Origin": "germline",
                }
            ],
            FIELDNAMES,
        )
        path = tmp_path / "clinvar.txt.gz"
        path.write_bytes(gz_bytes)
        data = load_clinvar_tsv(path)
        assert data["333"]["ClinicalSignificance"] == ""
        assert data["333"]["ReviewStatus"] == ""
        assert data["333"]["LastEvaluated"] == ""


# ---------------------------------------------------------------------------
# annotate_row
# ---------------------------------------------------------------------------


SAMPLE_CLINVAR_DATA = {
    "12345": {
        "ClinicalSignificance": "Pathogenic",
        "ReviewStatus": "reviewed by expert panel",
        "LastEvaluated": "2023-06-01",
    },
    "99999": {
        "ClinicalSignificance": "Benign",
        "ReviewStatus": "criteria provided, multiple submitters, no conflicts",
        "LastEvaluated": "2022-01-01",
    },
}


@pytest.mark.unit
class TestAnnotateRow:
    def _patch_resolve(self, mapping: dict[str, tuple[str, str]]):
        return patch(
            "annotate_clinvar.resolve_clinvar_ids",
            side_effect=lambda cid, cache, **kw: mapping.get(cid, ("", "")),
        )

    def test_empty_clingen_id(self):
        row = {"dna_clingen_allele_id": ""}
        out = annotate_row(row, SAMPLE_CLINVAR_DATA, {}, "clinvar.202601")
        assert out["clinvar.202601.clinical_significance"] == ""
        assert out["clinvar.202601.stars"] == ""

    def test_no_clingen_id_column(self):
        row = {}
        out = annotate_row(row, SAMPLE_CLINVAR_DATA, {}, "clinvar.202601")
        assert out["clinvar.202601.clinical_significance"] == ""

    def test_successful_annotation(self):
        row = {"dna_clingen_allele_id": "CA42"}
        with self._patch_resolve({"CA42": ("456", "12345")}):
            out = annotate_row(row, SAMPLE_CLINVAR_DATA, {}, "clinvar.202601")
        assert out["clinvar.202601.variation_id"] == "456"
        assert out["clinvar.202601.allele_id"] == "12345"
        assert out["clinvar.202601.clinical_significance"] == "Pathogenic"
        assert out["clinvar.202601.review_status"] == "reviewed by expert panel"
        assert out["clinvar.202601.stars"] == "3"
        assert out["clinvar.202601.last_review_date"] == "2023-06-01"

    def test_clingen_id_not_in_clinvar_tsv(self):
        row = {"dna_clingen_allele_id": "CA_MISS"}
        with self._patch_resolve({"CA_MISS": ("900", "00000")}):
            out = annotate_row(row, SAMPLE_CLINVAR_DATA, {}, "clinvar.202601")
        assert out["clinvar.202601.clinical_significance"] == ""
        assert out["clinvar.202601.variation_id"] == "900"
        assert out["clinvar.202601.allele_id"] == "00000"

    def test_clingen_id_resolves_to_empty(self):
        row = {"dna_clingen_allele_id": "CA_NONE"}
        with self._patch_resolve({"CA_NONE": ("", "")}):
            out = annotate_row(row, SAMPLE_CLINVAR_DATA, {}, "clinvar.202601")
        assert out["clinvar.202601.clinical_significance"] == ""
        assert out["clinvar.202601.variation_id"] == ""
        assert out["clinvar.202601.allele_id"] == ""

    def test_pipe_delimited_all_candidates_annotated(self):
        row = {"dna_clingen_allele_id": "CA_A|CA_B"}
        with self._patch_resolve({"CA_A": ("401", "12345"), "CA_B": ("402", "99999")}):
            out = annotate_row(row, SAMPLE_CLINVAR_DATA, {}, "clinvar.202601")
        assert out["clinvar.202601.variation_id"] == "401|402"
        assert out["clinvar.202601.allele_id"] == "12345|99999"
        assert out["clinvar.202601.clinical_significance"] == "Pathogenic|Benign"
        assert out["clinvar.202601.stars"] == "3|2"
        assert out["clinvar.202601.last_review_date"] == "2023-06-01|2022-01-01"

    def test_pipe_delimited_first_candidate_misses_second_hits(self):
        row = {"dna_clingen_allele_id": "CA_NOHIT|CA_B"}
        with self._patch_resolve({"CA_NOHIT": ("", ""), "CA_B": ("402", "99999")}):
            out = annotate_row(row, SAMPLE_CLINVAR_DATA, {}, "clinvar.202601")
        assert out["clinvar.202601.variation_id"] == "|402"
        assert out["clinvar.202601.allele_id"] == "|99999"
        assert out["clinvar.202601.clinical_significance"] == "|Benign"
        assert out["clinvar.202601.stars"] == "|2"

    def test_pipe_delimited_empty_slot_preserved(self):
        row = {"dna_clingen_allele_id": "CA_A||CA_B"}
        with self._patch_resolve({"CA_A": ("401", "12345"), "CA_B": ("402", "99999")}):
            out = annotate_row(row, SAMPLE_CLINVAR_DATA, {}, "clinvar.202601")
        assert out["clinvar.202601.variation_id"] == "401||402"
        assert out["clinvar.202601.allele_id"] == "12345||99999"
        assert out["clinvar.202601.clinical_significance"] == "Pathogenic||Benign"
        assert out["clinvar.202601.stars"] == "3||2"

    def test_custom_namespace_and_version(self):
        row = {"dna_clingen_allele_id": "CA42"}
        with self._patch_resolve({"CA42": ("456", "12345")}):
            out = annotate_row(row, SAMPLE_CLINVAR_DATA, {}, "cv.202501")
        assert "cv.202501.variation_id" in out
        assert out["cv.202501.variation_id"] == "456"
        assert "cv.202501.clinical_significance" in out
        assert out["cv.202501.clinical_significance"] == "Pathogenic"

    def test_custom_dna_column_name(self):
        row = {"my_dna_col": "CA42"}
        with self._patch_resolve({"CA42": ("456", "12345")}):
            out = annotate_row(
                row,
                SAMPLE_CLINVAR_DATA,
                {},
                "clinvar.202601",
                dna_clingen_allele_id_col="my_dna_col",
            )
        assert out["clinvar.202601.variation_id"] == "456"
        assert out["clinvar.202601.allele_id"] == "12345"
        assert out["clinvar.202601.clinical_significance"] == "Pathogenic"


@pytest.mark.integration
@pytest.mark.slow
def test_main_preserves_row_order_with_concurrency(tmp_path, monkeypatch):
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

    monkeypatch.setattr(mod, "fetch_clinvar_tsv", lambda y, m, c: Path("unused"))
    monkeypatch.setattr(mod, "load_clinvar_tsv", lambda p: {})

    def fake_annotate_row(
        row, clinvar_data, clingen_cache, col_prefix, dna_clingen_allele_id_col="dna_clingen_allele_id"
    ):
        import time as _time

        urn = row.get("variant_urn")
        if urn == "v1":
            _time.sleep(0.08)
        else:
            _time.sleep(0.02)
        return {
            f"{col_prefix}.variation_id": urn or "",
            f"{col_prefix}.allele_id": urn or "",
            f"{col_prefix}.clinical_significance": urn or "",
            f"{col_prefix}.review_status": "",
            f"{col_prefix}.stars": "",
            f"{col_prefix}.last_review_date": "",
        }

    monkeypatch.setattr(mod, "annotate_row", fake_annotate_row)

    mod.main([
        str(in_path),
        str(out_path),
        "--max-workers",
        "3",
        "--clinvar-version",
        "202601",
    ])

    with out_path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    assert [r["variant_urn"] for r in rows] == ["v1", "v2", "v3"]
    assert [r["clinvar.202601.clinical_significance"] for r in rows] == ["v1", "v2", "v3"]


@pytest.mark.integration
def test_main_invalid_max_workers_exits(tmp_path):
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

    with pytest.raises(SystemExit) as excinfo:
        mod.main([
            str(in_path),
            str(out_path),
            "--max-workers",
            "0",
        ])

    assert int(excinfo.value.code) == 1


@pytest.mark.integration
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

    monkeypatch.setattr(mod, "fetch_clinvar_tsv", lambda y, m, c: Path("unused"))
    monkeypatch.setattr(mod, "load_clinvar_tsv", lambda p: {})
    monkeypatch.setattr(
        mod,
        "annotate_row",
        lambda row, clinvar_data, clingen_cache, col_prefix, dna_clingen_allele_id_col="dna_clingen_allele_id": {
            f"{col_prefix}.variation_id": row["variant_urn"],
            f"{col_prefix}.allele_id": row["variant_urn"],
            f"{col_prefix}.clinical_significance": row["variant_urn"],
            f"{col_prefix}.review_status": "",
            f"{col_prefix}.stars": "",
            f"{col_prefix}.last_review_date": "",
        },
    )

    mod.main([
        str(in_path),
        str(out_path),
        "--clinvar-version",
        "202601",
        "--skip",
        "1",
        "--limit",
        "2",
    ])

    with out_path.open("r", encoding="utf-8", newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    header = out_path.read_text(encoding="utf-8").splitlines()[0].split("\t")
    id_cols = ["clinvar.202601.variation_id", "clinvar.202601.allele_id"]
    other_cols = [
        "clinvar.202601.clinical_significance",
        "clinvar.202601.review_status",
        "clinvar.202601.stars",
        "clinvar.202601.last_review_date",
    ]
    assert header.index(id_cols[0]) < header.index(other_cols[0])
    assert header.index(id_cols[1]) < header.index(other_cols[0])

    assert [r["variant_urn"] for r in rows] == ["v2", "v3"]
    assert [r["clinvar.202601.variation_id"] for r in rows] == ["v2", "v3"]
    assert [r["clinvar.202601.allele_id"] for r in rows] == ["v2", "v3"]
    assert [r["clinvar.202601.clinical_significance"] for r in rows] == ["v2", "v3"]
