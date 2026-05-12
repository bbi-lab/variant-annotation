"""Unit tests for src/annotate_missense_scores.py."""

from __future__ import annotations

import csv
from pathlib import Path
from unittest.mock import patch

import pytest

import src.annotate_missense_scores as mod
from src.annotate_missense_scores import (
    _snv_from_hgvs_g,
    _lookup_revel,
    _lookup_alphamissense,
    annotate_row,
    NC_TO_CHROM_GRCH38,
)


# ---------------------------------------------------------------------------
# _snv_from_hgvs_g
# ---------------------------------------------------------------------------

def test_snv_from_hgvs_g_recognises_snv():
    assert _snv_from_hgvs_g("NC_000001.11:g.69094T>A", NC_TO_CHROM_GRCH38) == ("1", 69094, "T", "A")


def test_snv_from_hgvs_g_lowercase_normalised():
    assert _snv_from_hgvs_g("NC_000007.14:g.117548628t>c", NC_TO_CHROM_GRCH38) == ("7", 117548628, "T", "C")


def test_snv_from_hgvs_g_rejects_deletion():
    assert _snv_from_hgvs_g("NC_000001.11:g.69094del", NC_TO_CHROM_GRCH38) is None


def test_snv_from_hgvs_g_rejects_mnv():
    # Two-base substitution is not a single SNV.
    assert _snv_from_hgvs_g("NC_000001.11:g.69094AT>GC", NC_TO_CHROM_GRCH38) is None


def test_snv_from_hgvs_g_rejects_unknown_accession():
    assert _snv_from_hgvs_g("NC_999999.99:g.100A>T", NC_TO_CHROM_GRCH38) is None


def test_snv_from_hgvs_g_empty_string():
    assert _snv_from_hgvs_g("", NC_TO_CHROM_GRCH38) is None


# ---------------------------------------------------------------------------
# _lookup_revel
# ---------------------------------------------------------------------------

_REVEL_LINE = "1\t69094\tT\tA\t0.4200"


def test_lookup_revel_found(tmp_path):
    """Monkeypatch _run_tabix to return a matching REVEL line."""
    dummy_path = tmp_path / "revel.tsv.gz"
    dummy_path.touch()
    cache: dict = {}

    with patch.object(mod, "_run_tabix", return_value=[_REVEL_LINE]):
        result = _lookup_revel(dummy_path, "1", 69094, "T", "A", cache)

    assert result == "0.4200"
    # Second call should use cache without calling _run_tabix again.
    with patch.object(mod, "_run_tabix", side_effect=AssertionError("should not call tabix")):
        assert _lookup_revel(dummy_path, "1", 69094, "T", "A", cache) == "0.4200"


def test_lookup_revel_not_found(tmp_path):
    dummy_path = tmp_path / "revel.tsv.gz"
    dummy_path.touch()
    cache: dict = {}

    with patch.object(mod, "_run_tabix", return_value=[]):
        result = _lookup_revel(dummy_path, "1", 69094, "T", "A", cache)

    assert result is None


def test_lookup_revel_takes_max_across_transcripts(tmp_path):
    """Multiple rows per position (different transcripts) → take max score."""
    lines = [
        "1\t69094\tT\tA\t0.3100",
        "1\t69094\tT\tA\t0.7800",  # higher — should win
        "1\t69094\tT\tA\t0.2000",
    ]
    dummy_path = tmp_path / "revel.tsv.gz"
    dummy_path.touch()
    cache: dict = {}

    with patch.object(mod, "_run_tabix", return_value=lines):
        result = _lookup_revel(dummy_path, "1", 69094, "T", "A", cache)

    assert result == "0.7800"


def test_lookup_revel_filters_ref_alt(tmp_path):
    """Lines that don't match ref/alt are ignored."""
    lines = [
        "1\t69094\tT\tC\t0.9900",  # different alt
        "1\t69094\tT\tA\t0.5500",  # correct
    ]
    dummy_path = tmp_path / "revel.tsv.gz"
    dummy_path.touch()
    cache: dict = {}

    with patch.object(mod, "_run_tabix", return_value=lines):
        result = _lookup_revel(dummy_path, "1", 69094, "T", "A", cache)

    assert result == "0.5500"


def test_lookup_revel_tries_chr_prefix(tmp_path):
    """If bare chrom returns nothing, tries 'chr' variant."""
    dummy_path = tmp_path / "revel.tsv.gz"
    dummy_path.touch()
    cache: dict = {}

    call_log: list[str] = []

    def fake_tabix(path: Path, chrom: str, pos: int) -> list[str]:
        call_log.append(chrom)
        if chrom == "chr1":
            return ["chr1\t69094\tT\tA\t0.3300"]
        return []

    with patch.object(mod, "_run_tabix", side_effect=fake_tabix):
        result = _lookup_revel(dummy_path, "1", 69094, "T", "A", cache)

    assert result == "0.3300"
    assert "1" in call_log
    assert "chr1" in call_log


# ---------------------------------------------------------------------------
# _lookup_alphamissense
# ---------------------------------------------------------------------------

_AM_LINE = (
    "chr1\t69094\tT\tA\thg38\tQ8NH21\tENST00000335137.4\tI1K\t0.7123\tlikely_pathogenic"
)


def test_lookup_alphamissense_found(tmp_path):
    dummy_path = tmp_path / "am.tsv.gz"
    dummy_path.touch()
    cache: dict = {}

    with patch.object(mod, "_run_tabix", return_value=[_AM_LINE]):
        result = _lookup_alphamissense(dummy_path, "chr1", 69094, "T", "A", cache)

    assert result == ("0.7123", "likely_pathogenic")


def test_lookup_alphamissense_not_found(tmp_path):
    dummy_path = tmp_path / "am.tsv.gz"
    dummy_path.touch()
    cache: dict = {}

    with patch.object(mod, "_run_tabix", return_value=[]):
        result = _lookup_alphamissense(dummy_path, "chr1", 69094, "T", "A", cache)

    assert result is None


def test_lookup_alphamissense_takes_max_pathogenicity(tmp_path):
    """Multiple transcript entries → return entry with highest pathogenicity."""
    lines = [
        "chr1\t69094\tT\tA\thg38\tP00001\tENST00000001\tI1K\t0.3000\tambiguous",
        "chr1\t69094\tT\tA\thg38\tP00001\tENST00000002\tI1K\t0.9100\tlikely_pathogenic",
        "chr1\t69094\tT\tA\thg38\tP00001\tENST00000003\tI1K\t0.1000\tlikely_benign",
    ]
    dummy_path = tmp_path / "am.tsv.gz"
    dummy_path.touch()
    cache: dict = {}

    with patch.object(mod, "_run_tabix", return_value=lines):
        result = _lookup_alphamissense(dummy_path, "chr1", 69094, "T", "A", cache)

    assert result == ("0.9100", "likely_pathogenic")


# ---------------------------------------------------------------------------
# annotate_row
# ---------------------------------------------------------------------------

def test_annotate_row_both_scores(tmp_path):
    revel_path = tmp_path / "revel.tsv.gz"
    am_path = tmp_path / "am.tsv.gz"
    revel_path.touch()
    am_path.touch()

    row = {"id": "v1", "mapped_hgvs_g": "NC_000001.11:g.69094T>A"}
    revel_cache: dict = {}
    am_cache: dict = {}

    def fake_tabix(path: Path, chrom: str, pos: int) -> list[str]:
        if "revel" in str(path):
            return ["1\t69094\tT\tA\t0.5500"]
        return ["chr1\t69094\tT\tA\thg38\tQ\tT\tI1K\t0.8200\tlikely_pathogenic"]

    with patch.object(mod, "_run_tabix", side_effect=fake_tabix):
        ann = annotate_row(
            row,
            nc_to_chrom=NC_TO_CHROM_GRCH38,
            mapped_hgvs_g_col="mapped_hgvs_g",
            revel_path=revel_path,
            alphamissense_path=am_path,
            revel_cache=revel_cache,
            am_cache=am_cache,
        )

    assert ann["revel.score"] == "0.5500"
    assert ann["alphamissense.pathogenicity"] == "0.8200"
    assert ann["alphamissense.class"] == "likely_pathogenic"


def test_annotate_row_non_snv_empty(tmp_path):
    """Indel HGVS → all scores empty."""
    revel_path = tmp_path / "revel.tsv.gz"
    revel_path.touch()

    row = {"mapped_hgvs_g": "NC_000001.11:g.69094del"}
    cache: dict = {}

    with patch.object(mod, "_run_tabix", side_effect=AssertionError("should not call tabix")):
        ann = annotate_row(
            row,
            nc_to_chrom=NC_TO_CHROM_GRCH38,
            mapped_hgvs_g_col="mapped_hgvs_g",
            revel_path=revel_path,
            alphamissense_path=None,
            revel_cache=cache,
            am_cache={},
        )

    assert ann["revel.score"] == ""


def test_annotate_row_pipe_aligned(tmp_path):
    """Pipe-delimited HGVS → pipe-delimited output aligned to positions."""
    revel_path = tmp_path / "revel.tsv.gz"
    revel_path.touch()

    row = {"mapped_hgvs_g": "NC_000001.11:g.69094T>A|NC_000001.11:g.69094T>C"}
    cache: dict = {}

    def fake_tabix(path: Path, chrom: str, pos: int) -> list[str]:
        # Return score only for T>A
        if chrom in ("1", "chr1"):
            return [
                "1\t69094\tT\tA\t0.6600",
                "1\t69094\tT\tC\t0.2200",
            ]
        return []

    with patch.object(mod, "_run_tabix", side_effect=fake_tabix):
        ann = annotate_row(
            row,
            nc_to_chrom=NC_TO_CHROM_GRCH38,
            mapped_hgvs_g_col="mapped_hgvs_g",
            revel_path=revel_path,
            alphamissense_path=None,
            revel_cache=cache,
            am_cache={},
        )

    parts = ann["revel.score"].split("|")
    assert len(parts) == 2
    assert parts[0] == "0.6600"
    assert parts[1] == "0.2200"


def test_annotate_row_empty_hgvs(tmp_path):
    """Empty mapped_hgvs_g column → empty scores, no tabix calls."""
    revel_path = tmp_path / "revel.tsv.gz"
    revel_path.touch()

    row = {"mapped_hgvs_g": ""}

    with patch.object(mod, "_run_tabix", side_effect=AssertionError("should not call tabix")):
        ann = annotate_row(
            row,
            nc_to_chrom=NC_TO_CHROM_GRCH38,
            mapped_hgvs_g_col="mapped_hgvs_g",
            revel_path=revel_path,
            alphamissense_path=None,
            revel_cache={},
            am_cache={},
        )

    assert ann["revel.score"] == ""


# ---------------------------------------------------------------------------
# main() integration
# ---------------------------------------------------------------------------

def _write_tsv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _read_tsv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def test_main_writes_score_columns(tmp_path, monkeypatch):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    revel_path = tmp_path / "revel.tsv.gz"
    am_path = tmp_path / "am.tsv.gz"
    revel_path.touch()
    am_path.touch()

    _write_tsv(
        in_path,
        [
            {"id": "v1", "mapped_hgvs_g": "NC_000001.11:g.69094T>A"},
            {"id": "v2", "mapped_hgvs_g": "NC_000001.11:g.69094del"},  # non-SNV → empty
        ],
        ["id", "mapped_hgvs_g"],
    )

    def fake_tabix(path: Path, chrom: str, pos: int) -> list[str]:
        if "revel" in str(path) and pos == 69094:
            return ["1\t69094\tT\tA\t0.7500"]
        if "am" in str(path) and pos == 69094:
            return ["chr1\t69094\tT\tA\thg38\tQ\tT\tI1K\t0.9000\tlikely_pathogenic"]
        return []

    monkeypatch.setattr(mod, "_run_tabix", fake_tabix)
    # Also stub out the tabix --version check.
    import subprocess
    monkeypatch.setattr(
        subprocess,
        "run",
        lambda *a, **kw: type("R", (), {"returncode": 0})()
        if "--version" in (a[0] if a else [])
        else subprocess.run.__wrapped__(*a, **kw)
        if hasattr(subprocess.run, "__wrapped__")
        else type("R", (), {"returncode": 0})(),
    )

    mod.main([
        str(in_path), str(out_path),
        "--revel-file", str(revel_path),
        "--alphamissense-file", str(am_path),
    ])

    rows = _read_tsv(out_path)
    assert len(rows) == 2
    assert rows[0]["revel.score"] == "0.7500"
    assert rows[0]["alphamissense.pathogenicity"] == "0.9000"
    assert rows[0]["alphamissense.class"] == "likely_pathogenic"
    # Non-SNV row → empty
    assert rows[1]["revel.score"] == ""
    assert rows[1]["alphamissense.pathogenicity"] == ""
