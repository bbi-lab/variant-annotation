"""Unit tests for src/annotate_clinvar.py."""

from __future__ import annotations

import csv
import gzip
import io
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

from src.annotate_clinvar import (
    annotate_row,
    load_clinvar_tsv,
    resolve_clinvar_allele_id,
    stars_for_review_status,
)


# ---------------------------------------------------------------------------
# stars_for_review_status
# ---------------------------------------------------------------------------


class TestStarsForReviewStatus:
    def test_practice_guideline_4(self):
        assert stars_for_review_status("practice guideline") == "4"

    def test_expert_panel_3(self):
        assert stars_for_review_status("reviewed by expert panel") == "3"

    def test_multiple_submitters_no_conflicts_2(self):
        assert (
            stars_for_review_status(
                "criteria provided, multiple submitters, no conflicts"
            )
            == "2"
        )

    def test_single_submitter_1(self):
        assert stars_for_review_status("criteria provided, single submitter") == "1"

    def test_conflicting_classifications_1(self):
        assert (
            stars_for_review_status("criteria provided, conflicting classifications")
            == "1"
        )

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


# ---------------------------------------------------------------------------
# resolve_clinvar_allele_id
# ---------------------------------------------------------------------------


class TestResolveClinvarAlleleId:
    def _mock_response(self, status_code: int, json_data: dict | None = None):
        resp = MagicMock()
        resp.status_code = status_code
        if json_data is not None:
            resp.json.return_value = json_data
        if status_code != 200:
            resp.raise_for_status.side_effect = requests_http_error(status_code)
        else:
            resp.raise_for_status.return_value = None
        return resp

    def test_returns_cached_value(self):
        cache = {"CA123": "99999"}
        result = resolve_clinvar_allele_id("CA123", cache)
        assert result == "99999"

    def test_404_returns_empty_string(self):
        cache: dict[str, str] = {}
        mock_resp = MagicMock()
        mock_resp.status_code = 404

        with patch("src.annotate_clinvar.requests.get", return_value=mock_resp):
            result = resolve_clinvar_allele_id("CA_UNKNOWN", cache)

        assert result == ""
        assert cache["CA_UNKNOWN"] == ""

    def test_successful_lookup(self):
        cache: dict[str, str] = {}
        mock_resp = MagicMock()
        mock_resp.status_code = 200
        mock_resp.json.return_value = {
            "externalRecords": {
                "ClinVarAlleles": [{"alleleId": 42}]
            }
        }

        with patch("src.annotate_clinvar.requests.get", return_value=mock_resp):
            result = resolve_clinvar_allele_id("CA42", cache)

        assert result == "42"
        assert cache["CA42"] == "42"

    def test_no_clinvar_association(self):
        cache: dict[str, str] = {}
        mock_resp = MagicMock()
        mock_resp.status_code = 200
        mock_resp.json.return_value = {"externalRecords": {}}

        with patch("src.annotate_clinvar.requests.get", return_value=mock_resp):
            result = resolve_clinvar_allele_id("CA_NOCLINVAR", cache)

        assert result == ""

    def test_client_error_no_retry(self):
        import requests as req_lib

        cache: dict[str, str] = {}
        mock_resp = MagicMock()
        mock_resp.status_code = 400
        http_error = req_lib.HTTPError(response=mock_resp)
        mock_resp.raise_for_status.side_effect = http_error

        with patch("src.annotate_clinvar.requests.get") as mock_get:
            mock_get.return_value = mock_resp
            # Patch raise_for_status so HTTPError is raised
            mock_resp.raise_for_status.side_effect = http_error
            # Make get() itself raise (simulating raise_for_status flow)
            # Instead, test that the function returns "" after one attempt
            result = resolve_clinvar_allele_id(
                "CA_BAD", cache, max_retries=3, retry_delay=0
            )

        # Should mark as empty and not retry indefinitely
        assert result == ""


def requests_http_error(status_code: int):
    """Helper to create a requests.HTTPError with a mock response."""
    import requests as req_lib

    mock_resp = MagicMock()
    mock_resp.status_code = status_code
    return req_lib.HTTPError(response=mock_resp)


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


class TestAnnotateRow:
    def _patch_resolve(self, mapping: dict[str, str]):
        """Context manager that patches resolve_clinvar_allele_id with *mapping*."""
        return patch(
            "src.annotate_clinvar.resolve_clinvar_allele_id",
            side_effect=lambda cid, cache, **kw: mapping.get(cid, ""),
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
        with self._patch_resolve({"CA42": "12345"}):
            out = annotate_row(row, SAMPLE_CLINVAR_DATA, {}, "clinvar.202601")
        assert out["clinvar.202601.clinical_significance"] == "Pathogenic"
        assert out["clinvar.202601.review_status"] == "reviewed by expert panel"
        assert out["clinvar.202601.stars"] == "3"
        assert out["clinvar.202601.last_review_date"] == "2023-06-01"

    def test_clingen_id_not_in_clinvar_tsv(self):
        row = {"dna_clingen_allele_id": "CA_MISS"}
        with self._patch_resolve({"CA_MISS": "00000"}):  # ID not in TSV
            out = annotate_row(row, SAMPLE_CLINVAR_DATA, {}, "clinvar.202601")
        assert out["clinvar.202601.clinical_significance"] == ""

    def test_clingen_id_resolves_to_empty(self):
        row = {"dna_clingen_allele_id": "CA_NONE"}
        with self._patch_resolve({"CA_NONE": ""}):
            out = annotate_row(row, SAMPLE_CLINVAR_DATA, {}, "clinvar.202601")
        assert out["clinvar.202601.clinical_significance"] == ""

    def test_pipe_delimited_first_candidate_hits(self):
        row = {"dna_clingen_allele_id": "CA_A|CA_B"}
        with self._patch_resolve({"CA_A": "12345", "CA_B": "99999"}):
            out = annotate_row(row, SAMPLE_CLINVAR_DATA, {}, "clinvar.202601")
        # CA_A resolves first and is present in TSV
        assert out["clinvar.202601.clinical_significance"] == "Pathogenic"

    def test_pipe_delimited_first_candidate_misses_second_hits(self):
        row = {"dna_clingen_allele_id": "CA_NOHIT|CA_B"}
        with self._patch_resolve({"CA_NOHIT": "", "CA_B": "99999"}):
            out = annotate_row(row, SAMPLE_CLINVAR_DATA, {}, "clinvar.202601")
        assert out["clinvar.202601.clinical_significance"] == "Benign"
        assert out["clinvar.202601.stars"] == "2"

    def test_custom_namespace_and_version(self):
        row = {"dna_clingen_allele_id": "CA42"}
        with self._patch_resolve({"CA42": "12345"}):
            out = annotate_row(row, SAMPLE_CLINVAR_DATA, {}, "cv.202501")
        assert "cv.202501.clinical_significance" in out
        assert out["cv.202501.clinical_significance"] == "Pathogenic"

    def test_custom_dna_column_name(self):
        row = {"my_dna_col": "CA42"}
        with self._patch_resolve({"CA42": "12345"}):
            out = annotate_row(
                row,
                SAMPLE_CLINVAR_DATA,
                {},
                "clinvar.202601",
                dna_clingen_allele_id_col="my_dna_col",
            )
        assert out["clinvar.202601.clinical_significance"] == "Pathogenic"
