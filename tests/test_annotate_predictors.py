"""Unit tests for src/annotate_predictors.py."""

from __future__ import annotations

import csv
from pathlib import Path
from unittest.mock import patch

import pytest

import gzip

import src.annotate_predictors as mod
from src.annotate_predictors import (
    _snv_from_hgvs_g,
    _lookup_revel,
    _lookup_alphamissense,
    _get_dbnsfp_col_indices,
    _lookup_mutpred2,
    _lookup_mutpred2_from_properties_file,
    _load_mutpred2_properties_file_cache,
    _lookup_revel_train,
    _lookup_mutpred2_train,
    _load_revel_training_variants,
    _load_mutpred2_training_variants,
    _unqualify_hgvs_p,
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
# _get_dbnsfp_col_indices
# ---------------------------------------------------------------------------

_DBNSFP_HEADER = (
    "#chr\tpos(1-based)\tref\talt\taaref\taaalt\trs_dbSNP\t"
    "MutPred2_score\tMutPred2_rankscore\tMutPred2_protID"
)


def test_get_dbnsfp_col_indices_parses_header(tmp_path):
    """Column names mapped to correct 0-based indices."""
    dummy = tmp_path / "db.gz"
    dummy.touch()

    # Clear module-level cache so the mock is actually called.
    mod._dbnsfp_col_index_cache.clear()

    fake_proc = type("P", (), {"stdout": _DBNSFP_HEADER + "\n", "returncode": 0})()
    with patch("subprocess.run", return_value=fake_proc):
        idx = _get_dbnsfp_col_indices(dummy)

    assert idx["chr"] == 0
    assert idx["pos(1-based)"] == 1
    assert idx["ref"] == 2
    assert idx["alt"] == 3
    assert idx["MutPred2_score"] == 7


def test_get_dbnsfp_col_indices_cached(tmp_path):
    """Second call returns cached value without calling subprocess."""
    dummy = tmp_path / "db2.gz"
    dummy.touch()
    mod._dbnsfp_col_index_cache.clear()

    fake_proc = type("P", (), {"stdout": _DBNSFP_HEADER + "\n", "returncode": 0})()
    with patch("subprocess.run", return_value=fake_proc) as mock_run:
        _get_dbnsfp_col_indices(dummy)
        _get_dbnsfp_col_indices(dummy)  # second call
        assert mock_run.call_count == 1  # subprocess called only once


def test_get_dbnsfp_col_indices_no_header_raises(tmp_path):
    dummy = tmp_path / "db3.gz"
    dummy.touch()
    mod._dbnsfp_col_index_cache.clear()

    fake_proc = type("P", (), {"stdout": "", "returncode": 0})()
    with patch("subprocess.run", return_value=fake_proc):
        with pytest.raises(ValueError, match="No header line"):
            _get_dbnsfp_col_indices(dummy)


# ---------------------------------------------------------------------------
# _lookup_mutpred2
# ---------------------------------------------------------------------------

_DBNSFP_IDX = {
    "chr": 0,
    "pos(1-based)": 1,
    "ref": 2,
    "alt": 3,
    "MutPred2_score": 4,
}


def _dbnsfp_line(chrom, pos, ref, alt, score):
    return f"{chrom}\t{pos}\t{ref}\t{alt}\t{score}"


def test_lookup_mutpred2_found(tmp_path):
    dummy = tmp_path / "db.gz"
    dummy.touch()
    cache: dict = {}
    mod._dbnsfp_col_index_cache[str(dummy)] = _DBNSFP_IDX

    with patch.object(mod, "_run_tabix", return_value=[_dbnsfp_line("1", 69094, "T", "A", "0.6700")]):
        result = _lookup_mutpred2(dummy, "1", 69094, "T", "A", cache)

    assert result == "0.6700"


def test_lookup_mutpred2_not_found(tmp_path):
    dummy = tmp_path / "db.gz"
    dummy.touch()
    cache: dict = {}
    mod._dbnsfp_col_index_cache[str(dummy)] = _DBNSFP_IDX

    with patch.object(mod, "_run_tabix", return_value=[]):
        result = _lookup_mutpred2(dummy, "1", 69094, "T", "A", cache)

    assert result is None


def test_lookup_mutpred2_semicolon_max(tmp_path):
    """Multiple transcript scores (semicolon-separated) → take max."""
    dummy = tmp_path / "db.gz"
    dummy.touch()
    cache: dict = {}
    mod._dbnsfp_col_index_cache[str(dummy)] = _DBNSFP_IDX

    line = _dbnsfp_line("1", 69094, "T", "A", "0.3100;0.8900;0.1200")
    with patch.object(mod, "_run_tabix", return_value=[line]):
        result = _lookup_mutpred2(dummy, "1", 69094, "T", "A", cache)

    assert result == "0.8900"


def test_lookup_mutpred2_dot_is_missing(tmp_path):
    """Periods and 'NA' are treated as no score."""
    dummy = tmp_path / "db.gz"
    dummy.touch()
    cache: dict = {}
    mod._dbnsfp_col_index_cache[str(dummy)] = _DBNSFP_IDX

    line = _dbnsfp_line("1", 69094, "T", "A", ".")
    with patch.object(mod, "_run_tabix", return_value=[line]):
        result = _lookup_mutpred2(dummy, "1", 69094, "T", "A", cache)

    assert result is None


def test_lookup_mutpred2_uses_cache(tmp_path):
    dummy = tmp_path / "db.gz"
    dummy.touch()
    cache: dict = {}
    mod._dbnsfp_col_index_cache[str(dummy)] = _DBNSFP_IDX

    with patch.object(mod, "_run_tabix", return_value=[_dbnsfp_line("1", 69094, "T", "A", "0.5000")]):
        _lookup_mutpred2(dummy, "1", 69094, "T", "A", cache)

    with patch.object(mod, "_run_tabix", side_effect=AssertionError("should not call tabix")):
        result = _lookup_mutpred2(dummy, "1", 69094, "T", "A", cache)

    assert result == "0.5000"


# ---------------------------------------------------------------------------
# _load_mutpred2_properties_file_cache / _lookup_mutpred2_from_properties_file
# ---------------------------------------------------------------------------

_MP2_PROPERTIES_HEADER = (
    "Dataset,Gene,mavedb_variant_urn,Chrom,Strand,hg38_start,hg38_end,"
    "ref_allele,alt_allele,AA,protein_id,transcript_id,gene_id,gene_symbol,"
    "Substitution,MutPred2 score,Mechanisms"
)


def _mp2_line(urn, chrom, start, stop, ref, alt, score):
    return (
        f"ds,GENE,{urn},{chrom},-1.0,{start},{stop},{ref},{alt},T2A,"
        f"ENSP1,ENST1,ENSG1,GENE,T2A,{score},[]"
    )


def _write_mp2_properties_file(path: Path, lines: list[str], gzipped: bool) -> Path:
    text = _MP2_PROPERTIES_HEADER + "\n" + "\n".join(lines) + "\n"
    if gzipped:
        file_path = path / "mp2.csv.gz"
        with gzip.open(file_path, "wt", encoding="utf-8") as fh:
            fh.write(text)
    else:
        file_path = path / "mp2.csv"
        file_path.write_text(text, encoding="utf-8")
    return file_path


def test_load_mutpred2_properties_file_cache_plain_csv(tmp_path):
    lines = [_mp2_line("urn:mavedb:1#1", "17", 100, 100, "A", "G", "0.5000")]
    file_path = _write_mp2_properties_file(tmp_path, lines, gzipped=False)

    cache = _load_mutpred2_properties_file_cache(str(file_path))

    assert cache[("urn:mavedb:1#1", "17", 100, 100, "A", "G")] == "0.5000"


def test_load_mutpred2_properties_file_cache_gzipped(tmp_path):
    lines = [_mp2_line("urn:mavedb:1#1", "17", 100, 100, "A", "G", "0.5000")]
    file_path = _write_mp2_properties_file(tmp_path, lines, gzipped=True)

    cache = _load_mutpred2_properties_file_cache(str(file_path))

    assert cache[("urn:mavedb:1#1", "17", 100, 100, "A", "G")] == "0.5000"


def test_load_mutpred2_properties_file_cache_multiple_candidates_same_urn(tmp_path):
    """One variant_urn can repeat across multiple DNA reverse-translation rows."""
    lines = [
        _mp2_line("urn:mavedb:1#1", "17", 100, 102, "ACT", "GCA", "0.3000"),
        _mp2_line("urn:mavedb:1#1", "17", 100, 102, "ACT", "GCC", "0.3000"),
    ]
    file_path = _write_mp2_properties_file(tmp_path, lines, gzipped=False)

    cache = _load_mutpred2_properties_file_cache(str(file_path))

    assert cache[("urn:mavedb:1#1", "17", 100, 102, "ACT", "GCA")] == "0.3000"
    assert cache[("urn:mavedb:1#1", "17", 100, 102, "ACT", "GCC")] == "0.3000"


def test_load_mutpred2_properties_file_cache_missing_column_raises(tmp_path):
    file_path = tmp_path / "bad.csv"
    file_path.write_text("Dataset,Gene,mavedb_variant_urn\nds,GENE,urn:1\n", encoding="utf-8")

    with pytest.raises(ValueError, match="missing column"):
        _load_mutpred2_properties_file_cache(str(file_path))


def test_load_mutpred2_properties_file_cache_skips_bad_rows(tmp_path):
    lines = [
        _mp2_line("urn:mavedb:1#1", "17", 100, 100, "A", "G", "0.5000"),
        _mp2_line("", "17", 100, 100, "A", "G", "0.5000"),  # blank urn
        _mp2_line("urn:mavedb:1#2", "17", "notanumber", 100, "A", "G", "0.5000"),  # bad start
        _mp2_line("urn:mavedb:1#3", "17", 100, 100, "A", "G", ""),  # blank score
    ]
    file_path = _write_mp2_properties_file(tmp_path, lines, gzipped=False)

    cache = _load_mutpred2_properties_file_cache(str(file_path))

    assert len(cache) == 1
    assert cache[("urn:mavedb:1#1", "17", 100, 100, "A", "G")] == "0.5000"


def test_lookup_mutpred2_from_properties_file_found():
    cache = {("urn:mavedb:1#1", "17", 100, 100, "A", "G"): "0.5000"}
    result = _lookup_mutpred2_from_properties_file(cache, "urn:mavedb:1#1", "17", "100", "100", "A", "G")
    assert result == "0.5000"


def test_lookup_mutpred2_from_properties_file_not_found():
    cache = {("urn:mavedb:1#1", "17", 100, 100, "A", "G"): "0.5000"}
    result = _lookup_mutpred2_from_properties_file(cache, "urn:mavedb:1#1", "17", "100", "100", "A", "C")
    assert result is None


def test_lookup_mutpred2_from_properties_file_tries_chr_prefix():
    cache = {("urn:mavedb:1#1", "chr17", 100, 100, "A", "G"): "0.5000"}
    result = _lookup_mutpred2_from_properties_file(cache, "urn:mavedb:1#1", "17", "100", "100", "A", "G")
    assert result == "0.5000"


def test_lookup_mutpred2_from_properties_file_missing_field_returns_none():
    cache = {("urn:mavedb:1#1", "17", 100, 100, "A", "G"): "0.5000"}
    result = _lookup_mutpred2_from_properties_file(cache, "", "17", "100", "100", "A", "G")
    assert result is None


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
    dbnsfp_path = tmp_path / "db.tsv.gz"
    revel_path.touch()
    am_path.touch()
    dbnsfp_path.touch()
    mod._dbnsfp_col_index_cache[str(dbnsfp_path)] = _DBNSFP_IDX

    row = {"id": "v1", "mapped_hgvs_g": "NC_000001.11:g.69094T>A"}
    revel_cache: dict = {}
    am_cache: dict = {}
    mutpred2_cache: dict = {}

    def fake_tabix(path: Path, chrom: str, pos: int) -> list[str]:
        if "revel" in str(path):
            return ["1\t69094\tT\tA\t0.5500"]
        if "db" in str(path):
            return [_dbnsfp_line("1", 69094, "T", "A", "0.7200")]
        return ["chr1\t69094\tT\tA\thg38\tQ\tT\tI1K\t0.8200\tlikely_pathogenic"]

    with patch.object(mod, "_run_tabix", side_effect=fake_tabix):
        ann = annotate_row(
            row,
            nc_to_chrom=NC_TO_CHROM_GRCH38,
            mapped_hgvs_g_col="mapped_hgvs_g",
            revel_path=revel_path,
            alphamissense_path=am_path,
            dbnsfp_path=dbnsfp_path,
            revel_cache=revel_cache,
            am_cache=am_cache,
            mutpred2_cache=mutpred2_cache,
        )

    assert ann["revel.score"] == "0.5500"
    assert ann["alphamissense.pathogenicity"] == "0.8200"
    assert ann["alphamissense.class"] == "likely_pathogenic"
    assert ann["mutpred2.score"] == "0.7200"


def test_annotate_row_with_dbnsfp_only(tmp_path):
    """dbnsfp_path without revel/am → only mutpred2.score in output (single value)."""
    dbnsfp_path = tmp_path / "db.tsv.gz"
    dbnsfp_path.touch()
    mod._dbnsfp_col_index_cache[str(dbnsfp_path)] = _DBNSFP_IDX

    row = {"mapped_hgvs_g": "NC_000001.11:g.69094T>A"}

    with patch.object(mod, "_run_tabix", return_value=[_dbnsfp_line("1", 69094, "T", "A", "0.6300")]):
        ann = annotate_row(
            row,
            nc_to_chrom=NC_TO_CHROM_GRCH38,
            mapped_hgvs_g_col="mapped_hgvs_g",
            revel_path=None,
            alphamissense_path=None,
            dbnsfp_path=dbnsfp_path,
            revel_cache={},
            am_cache={},
        )

    assert "mutpred2.score" in ann
    assert ann["mutpred2.score"] == "0.6300"
    assert "revel.score" not in ann


def test_annotate_row_mutpred2_single_score_across_candidates(tmp_path):
    """Pipe-delimited candidates → MutPred2 emits one score (protein model).

    Candidates are different DNA spellings of the same amino acid change, so
    mutpred2.score must be a plain string, not pipe-delimited.
    """
    revel_path = tmp_path / "revel.tsv.gz"
    dbnsfp_path = tmp_path / "db.tsv.gz"
    revel_path.touch()
    dbnsfp_path.touch()
    mod._dbnsfp_col_index_cache[str(dbnsfp_path)] = _DBNSFP_IDX

    # Two reverse-translation candidates for the same protein change.
    row = {"mapped_hgvs_g": "NC_000001.11:g.69094T>A|NC_000001.11:g.69094T>C"}

    def fake_tabix(path: Path, chrom: str, pos: int) -> list[str]:
        if "revel" in str(path):
            return [
                "1\t69094\tT\tA\t0.5500",
                "1\t69094\tT\tC\t0.3300",
            ]
        # dbNSFP: first candidate has score 0.6700, second has 0.8100.
        return [
            _dbnsfp_line("1", 69094, "T", "A", "0.6700"),
            _dbnsfp_line("1", 69094, "T", "C", "0.8100"),
        ]

    with patch.object(mod, "_run_tabix", side_effect=fake_tabix):
        ann = annotate_row(
            row,
            nc_to_chrom=NC_TO_CHROM_GRCH38,
            mapped_hgvs_g_col="mapped_hgvs_g",
            revel_path=revel_path,
            alphamissense_path=None,
            dbnsfp_path=dbnsfp_path,
            revel_cache={},
            am_cache={},
        )

    # REVEL is DNA-level → pipe-aligned.
    assert ann["revel.score"] == "0.5500|0.3300"
    # MutPred2 is protein-level → single best score, no pipe.
    assert "|" not in ann["mutpred2.score"]
    assert ann["mutpred2.score"] == "0.8100"


def test_annotate_row_mutpred2_properties_file_pipe_aligned_per_candidate(tmp_path):
    """Alt source: MutPred2 is looked up per DNA RT candidate, pipe-aligned."""
    mutpred2_properties_cache = {
        ("urn:mavedb:1#1", "17", 100, 100, "A", "G"): "0.4000",
        ("urn:mavedb:1#1", "17", 100, 100, "A", "T"): "0.9000",
    }
    row = {
        "variant_urn": "urn:mavedb:1#1",
        "mapped_hgvs_g": "NC_000017.11:g.100A>G|NC_000017.11:g.100A>T",
        "mapped_hgvs_g_chromosome": "17|17",
        "mapped_hgvs_g_start": "100|100",
        "mapped_hgvs_g_stop": "100|100",
        "mapped_hgvs_g_ref": "A|A",
        "mapped_hgvs_g_alt": "G|T",
    }

    ann = annotate_row(
        row,
        nc_to_chrom=NC_TO_CHROM_GRCH38,
        mapped_hgvs_g_col="mapped_hgvs_g",
        revel_path=None,
        alphamissense_path=None,
        dbnsfp_path=None,
        revel_cache={},
        am_cache={},
        mutpred2_properties_cache=mutpred2_properties_cache,
        variant_urn_col="variant_urn",
        mapped_hgvs_g_chromosome_col="mapped_hgvs_g_chromosome",
        mapped_hgvs_g_start_col="mapped_hgvs_g_start",
        mapped_hgvs_g_stop_col="mapped_hgvs_g_stop",
        mapped_hgvs_g_ref_col="mapped_hgvs_g_ref",
        mapped_hgvs_g_alt_col="mapped_hgvs_g_alt",
    )

    assert ann["mutpred2.score"] == "0.4000|0.9000"


def test_annotate_row_mutpred2_properties_file_missing_candidate_is_empty(tmp_path):
    """A candidate absent from the properties file produces an empty slot, not a dropped one."""
    mutpred2_properties_cache = {
        ("urn:mavedb:1#1", "17", 100, 100, "A", "G"): "0.4000",
    }
    row = {
        "variant_urn": "urn:mavedb:1#1",
        "mapped_hgvs_g": "NC_000017.11:g.100A>G|NC_000017.11:g.100A>T",
        "mapped_hgvs_g_chromosome": "17|17",
        "mapped_hgvs_g_start": "100|100",
        "mapped_hgvs_g_stop": "100|100",
        "mapped_hgvs_g_ref": "A|A",
        "mapped_hgvs_g_alt": "G|T",
    }

    ann = annotate_row(
        row,
        nc_to_chrom=NC_TO_CHROM_GRCH38,
        mapped_hgvs_g_col="mapped_hgvs_g",
        revel_path=None,
        alphamissense_path=None,
        dbnsfp_path=None,
        revel_cache={},
        am_cache={},
        mutpred2_properties_cache=mutpred2_properties_cache,
        variant_urn_col="variant_urn",
        mapped_hgvs_g_chromosome_col="mapped_hgvs_g_chromosome",
        mapped_hgvs_g_start_col="mapped_hgvs_g_start",
        mapped_hgvs_g_stop_col="mapped_hgvs_g_stop",
        mapped_hgvs_g_ref_col="mapped_hgvs_g_ref",
        mapped_hgvs_g_alt_col="mapped_hgvs_g_alt",
    )

    assert ann["mutpred2.score"] == "0.4000|"


def test_annotate_row_mutpred2_properties_file_takes_precedence_over_dbnsfp(tmp_path):
    """When both sources are configured, the properties file wins and dbNSFP is skipped."""
    dbnsfp_path = tmp_path / "db.tsv.gz"
    dbnsfp_path.touch()
    mod._dbnsfp_col_index_cache[str(dbnsfp_path)] = _DBNSFP_IDX

    mutpred2_properties_cache = {
        ("urn:mavedb:1#1", "17", 100, 100, "A", "G"): "0.4000",
    }
    row = {
        "variant_urn": "urn:mavedb:1#1",
        "mapped_hgvs_g": "NC_000017.11:g.100A>G",
        "mapped_hgvs_g_chromosome": "17",
        "mapped_hgvs_g_start": "100",
        "mapped_hgvs_g_stop": "100",
        "mapped_hgvs_g_ref": "A",
        "mapped_hgvs_g_alt": "G",
    }

    with patch.object(mod, "_run_tabix", side_effect=AssertionError("dbNSFP should not be queried")):
        ann = annotate_row(
            row,
            nc_to_chrom=NC_TO_CHROM_GRCH38,
            mapped_hgvs_g_col="mapped_hgvs_g",
            revel_path=None,
            alphamissense_path=None,
            dbnsfp_path=dbnsfp_path,
            revel_cache={},
            am_cache={},
            mutpred2_properties_cache=mutpred2_properties_cache,
            variant_urn_col="variant_urn",
            mapped_hgvs_g_chromosome_col="mapped_hgvs_g_chromosome",
            mapped_hgvs_g_start_col="mapped_hgvs_g_start",
            mapped_hgvs_g_stop_col="mapped_hgvs_g_stop",
            mapped_hgvs_g_ref_col="mapped_hgvs_g_ref",
            mapped_hgvs_g_alt_col="mapped_hgvs_g_alt",
        )

    assert ann["mutpred2.score"] == "0.4000"


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


def test_main_mutpred2_properties_file_pipe_aligned_output(tmp_path):
    """End-to-end: protein-level row with multiple RT candidates → pipe-aligned mutpred2.score."""
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    mp2_path = tmp_path / "mp2.csv"
    mp2_path.write_text(
        _MP2_PROPERTIES_HEADER + "\n"
        + _mp2_line("urn:mavedb:1#1", "17", 100, 100, "A", "G", "0.4000") + "\n"
        + _mp2_line("urn:mavedb:1#1", "17", 100, 100, "A", "T", "0.9000") + "\n",
        encoding="utf-8",
    )

    _write_tsv(
        in_path,
        [
            {
                "variant_urn": "urn:mavedb:1#1",
                "mapped_hgvs_g": "NC_000017.11:g.100A>G|NC_000017.11:g.100A>T",
                "mapped_hgvs_g_chromosome": "17|17",
                "mapped_hgvs_g_start": "100|100",
                "mapped_hgvs_g_stop": "100|100",
                "mapped_hgvs_g_ref": "A|A",
                "mapped_hgvs_g_alt": "G|T",
            },
        ],
        [
            "variant_urn", "mapped_hgvs_g", "mapped_hgvs_g_chromosome",
            "mapped_hgvs_g_start", "mapped_hgvs_g_stop", "mapped_hgvs_g_ref", "mapped_hgvs_g_alt",
        ],
    )

    mod.main([
        str(in_path), str(out_path),
        "--mutpred2-properties-file", str(mp2_path),
    ])

    rows = _read_tsv(out_path)
    assert len(rows) == 1
    assert rows[0]["mutpred2.score"] == "0.4000|0.9000"


def test_main_mutpred2_properties_file_takes_precedence_over_dbnsfp(tmp_path, monkeypatch):
    """--dbnsfp-file is ignored for mutpred2.score when --mutpred2-properties-file is also given."""
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    mp2_path = tmp_path / "mp2.csv"
    dbnsfp_path = tmp_path / "db.tsv.gz"
    dbnsfp_path.touch()
    mp2_path.write_text(
        _MP2_PROPERTIES_HEADER + "\n"
        + _mp2_line("urn:mavedb:1#1", "17", 100, 100, "A", "G", "0.4000") + "\n",
        encoding="utf-8",
    )

    _write_tsv(
        in_path,
        [
            {
                "variant_urn": "urn:mavedb:1#1",
                "mapped_hgvs_g": "NC_000017.11:g.100A>G",
                "mapped_hgvs_g_chromosome": "17",
                "mapped_hgvs_g_start": "100",
                "mapped_hgvs_g_stop": "100",
                "mapped_hgvs_g_ref": "A",
                "mapped_hgvs_g_alt": "G",
            },
        ],
        [
            "variant_urn", "mapped_hgvs_g", "mapped_hgvs_g_chromosome",
            "mapped_hgvs_g_start", "mapped_hgvs_g_stop", "mapped_hgvs_g_ref", "mapped_hgvs_g_alt",
        ],
    )

    monkeypatch.setattr(mod, "_run_tabix", lambda *a, **kw: (_ for _ in ()).throw(
        AssertionError("dbNSFP should not be queried")
    ))
    import subprocess
    monkeypatch.setattr(subprocess, "run", lambda *a, **kw: type("R", (), {"returncode": 0})())

    mod.main([
        str(in_path), str(out_path),
        "--mutpred2-properties-file", str(mp2_path),
        "--dbnsfp-file", str(dbnsfp_path),
    ])

    rows = _read_tsv(out_path)
    assert rows[0]["mutpred2.score"] == "0.4000"


# ---------------------------------------------------------------------------
# _unqualify_hgvs_p
# ---------------------------------------------------------------------------

def test_unqualify_hgvs_p_strips_accession():
    assert _unqualify_hgvs_p("NP_000546.1:p.Asn1643His") == "p.Asn1643His"


def test_unqualify_hgvs_p_passthrough_when_no_colon():
    assert _unqualify_hgvs_p("p.Asn1643His") == "p.Asn1643His"


def test_unqualify_hgvs_p_strips_whitespace():
    assert _unqualify_hgvs_p("  NP_000546.1:p.Asn1643His  ") == "p.Asn1643His"


# ---------------------------------------------------------------------------
# _load_revel_training_variants
# ---------------------------------------------------------------------------

def _write_revel_training_file(path: Path, rows: list[dict]) -> Path:
    file_path = path / "revel_training.tsv"
    fieldnames = ["gene_symbol", "unqualified_hgvs_p", "chromosome", "hg38_start", "hg38_end", "ref_allele", "alt_allele"]
    _write_tsv(file_path, rows, fieldnames)
    return file_path


def test_load_revel_training_variants_basic(tmp_path):
    file_path = _write_revel_training_file(tmp_path, [
        {
            "gene_symbol": "ASPA", "unqualified_hgvs_p": "p.Cys4Arg",
            "chromosome": "17", "hg38_start": "3476169.0", "hg38_end": "3476169.0",
            "ref_allele": "T", "alt_allele": "C",
        },
    ])

    training_set = _load_revel_training_variants(str(file_path))

    assert training_set == {("17", 3476169, 3476169, "T", "C")}


def test_load_revel_training_variants_missing_column_raises(tmp_path):
    file_path = tmp_path / "bad.tsv"
    file_path.write_text("gene_symbol\tunqualified_hgvs_p\nASPA\tp.Cys4Arg\n", encoding="utf-8")

    with pytest.raises(ValueError, match="missing column"):
        _load_revel_training_variants(str(file_path))


def test_load_revel_training_variants_skips_bad_positions(tmp_path):
    file_path = _write_revel_training_file(tmp_path, [
        {
            "gene_symbol": "ASPA", "unqualified_hgvs_p": "p.Cys4Arg",
            "chromosome": "17", "hg38_start": "notanumber", "hg38_end": "3476169.0",
            "ref_allele": "T", "alt_allele": "C",
        },
    ])

    training_set = _load_revel_training_variants(str(file_path))

    assert training_set == set()


# ---------------------------------------------------------------------------
# _load_mutpred2_training_variants
# ---------------------------------------------------------------------------

def test_load_mutpred2_training_variants_unqualified_schema(tmp_path):
    file_path = tmp_path / "mp2_training.tsv"
    _write_tsv(
        file_path,
        [{"gene_symbol": "TSC2", "unqualified_hgvs_p": "p.Asn1643His"}],
        ["gene_symbol", "unqualified_hgvs_p"],
    )

    schema, training_set = _load_mutpred2_training_variants(str(file_path))

    assert schema == "unqualified"
    assert training_set == {("TSC2", "p.Asn1643His")}


def test_load_mutpred2_training_variants_qualified_schema(tmp_path):
    file_path = tmp_path / "mp2_training.tsv"
    _write_tsv(
        file_path,
        [{"hgvs_p": "NP_000546.1:p.Asn1643His"}],
        ["hgvs_p"],
    )

    schema, training_set = _load_mutpred2_training_variants(str(file_path))

    assert schema == "qualified"
    assert training_set == {"NP_000546.1:p.Asn1643His"}


def test_load_mutpred2_training_variants_missing_columns_raises(tmp_path):
    file_path = tmp_path / "bad.tsv"
    file_path.write_text("gene_symbol\nTSC2\n", encoding="utf-8")

    with pytest.raises(ValueError, match="must have either"):
        _load_mutpred2_training_variants(str(file_path))


# ---------------------------------------------------------------------------
# _lookup_revel_train
# ---------------------------------------------------------------------------

def test_lookup_revel_train_found():
    training_set = {("17", 3476169, 3476169, "T", "C")}
    assert _lookup_revel_train(training_set, "17", "3476169", "3476169", "T", "C") is True


def test_lookup_revel_train_not_found():
    training_set = {("17", 3476169, 3476169, "T", "C")}
    assert _lookup_revel_train(training_set, "17", "100", "100", "A", "G") is False


def test_lookup_revel_train_tries_chr_prefix():
    training_set = {("17", 3476169, 3476169, "T", "C")}
    assert _lookup_revel_train(training_set, "chr17", "3476169", "3476169", "T", "C") is True


def test_lookup_revel_train_missing_field_is_false():
    training_set = {("17", 3476169, 3476169, "T", "C")}
    assert _lookup_revel_train(training_set, "17", "", "3476169", "T", "C") is False


# ---------------------------------------------------------------------------
# _lookup_mutpred2_train
# ---------------------------------------------------------------------------

def test_lookup_mutpred2_train_unqualified_match():
    training_set = {("TSC2", "p.Asn1643His")}
    assert _lookup_mutpred2_train("unqualified", training_set, "TSC2", "NP_001071.3:p.Asn1643His") is True


def test_lookup_mutpred2_train_unqualified_no_gene_match():
    training_set = {("TSC2", "p.Asn1643His")}
    assert _lookup_mutpred2_train("unqualified", training_set, "OTHERGENE", "p.Asn1643His") is False


def test_lookup_mutpred2_train_qualified_match_ignores_gene():
    training_set = {"NP_001071.3:p.Asn1643His"}
    assert _lookup_mutpred2_train("qualified", training_set, "", "NP_001071.3:p.Asn1643His") is True


def test_lookup_mutpred2_train_pipe_delimited_any_segment_matches():
    training_set = {("TSC2", "p.Asn1643His")}
    mapped_hgvs_p = "NP_999999.1:p.Gln2Ter|NP_001071.3:p.Asn1643His"
    assert _lookup_mutpred2_train("unqualified", training_set, "TSC2", mapped_hgvs_p) is True


def test_lookup_mutpred2_train_empty_value_is_false():
    training_set = {("TSC2", "p.Asn1643His")}
    assert _lookup_mutpred2_train("unqualified", training_set, "TSC2", "") is False


# ---------------------------------------------------------------------------
# annotate_row: revel.train / mutpred2.train
# ---------------------------------------------------------------------------

def test_annotate_row_revel_train_pipe_aligned_per_candidate(tmp_path):
    revel_training_set = {("17", 100, 100, "A", "G")}
    row = {
        "mapped_hgvs_g": "NC_000017.11:g.100A>G|NC_000017.11:g.100A>T",
        "mapped_hgvs_g_chromosome": "17|17",
        "mapped_hgvs_g_start": "100|100",
        "mapped_hgvs_g_stop": "100|100",
        "mapped_hgvs_g_ref": "A|A",
        "mapped_hgvs_g_alt": "G|T",
    }

    ann = annotate_row(
        row,
        nc_to_chrom=NC_TO_CHROM_GRCH38,
        mapped_hgvs_g_col="mapped_hgvs_g",
        revel_path=None,
        alphamissense_path=None,
        revel_cache={},
        am_cache={},
        mapped_hgvs_g_chromosome_col="mapped_hgvs_g_chromosome",
        mapped_hgvs_g_start_col="mapped_hgvs_g_start",
        mapped_hgvs_g_stop_col="mapped_hgvs_g_stop",
        mapped_hgvs_g_ref_col="mapped_hgvs_g_ref",
        mapped_hgvs_g_alt_col="mapped_hgvs_g_alt",
        revel_training_set=revel_training_set,
    )

    assert ann["revel.train"] == "true|false"


def test_annotate_row_mutpred2_train_duplicated_across_candidates(tmp_path):
    """MutPred2 training is protein-level: same true/false repeats for every DNA candidate."""
    mutpred2_training_set = {("TSC2", "p.Asn1643His")}
    row = {
        "mapped_hgvs_g": "NC_000016.10:g.100A>G|NC_000016.10:g.200A>T",
        "gene_symbol": "TSC2",
        "mapped_hgvs_p": "NP_001071.3:p.Asn1643His",
    }

    ann = annotate_row(
        row,
        nc_to_chrom=NC_TO_CHROM_GRCH38,
        mapped_hgvs_g_col="mapped_hgvs_g",
        revel_path=None,
        alphamissense_path=None,
        revel_cache={},
        am_cache={},
        mutpred2_training_schema="unqualified",
        mutpred2_training_set=mutpred2_training_set,
        gene_symbol_col="gene_symbol",
        mapped_hgvs_p_col="mapped_hgvs_p",
    )

    assert ann["mutpred2.train"] == "true|true"


def test_annotate_row_mutpred2_train_single_candidate_no_pipe(tmp_path):
    mutpred2_training_set = {("TSC2", "p.Asn1643His")}
    row = {
        "mapped_hgvs_g": "NC_000016.10:g.100A>G",
        "gene_symbol": "OTHERGENE",
        "mapped_hgvs_p": "NP_001071.3:p.Asn1643His",
    }

    ann = annotate_row(
        row,
        nc_to_chrom=NC_TO_CHROM_GRCH38,
        mapped_hgvs_g_col="mapped_hgvs_g",
        revel_path=None,
        alphamissense_path=None,
        revel_cache={},
        am_cache={},
        mutpred2_training_schema="unqualified",
        mutpred2_training_set=mutpred2_training_set,
        gene_symbol_col="gene_symbol",
        mapped_hgvs_p_col="mapped_hgvs_p",
    )

    assert ann["mutpred2.train"] == "false"


# ---------------------------------------------------------------------------
# main(): revel.train / mutpred2.train end-to-end
# ---------------------------------------------------------------------------

def test_main_training_files_only_writes_train_columns(tmp_path):
    """--revel-training-file and --mutpred2-training-file work without any score source."""
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"

    revel_training_path = _write_revel_training_file(tmp_path, [
        {
            "gene_symbol": "ASPA", "unqualified_hgvs_p": "p.Cys4Arg",
            "chromosome": "17", "hg38_start": "100", "hg38_end": "100",
            "ref_allele": "A", "alt_allele": "G",
        },
    ])
    mutpred2_training_path = tmp_path / "mp2_training.tsv"
    _write_tsv(
        mutpred2_training_path,
        [{"gene_symbol": "TSC2", "unqualified_hgvs_p": "p.Asn1643His"}],
        ["gene_symbol", "unqualified_hgvs_p"],
    )

    _write_tsv(
        in_path,
        [
            {
                "gene_symbol": "ASPA",
                "mapped_hgvs_g": "NC_000017.11:g.100A>G",
                "mapped_hgvs_g_chromosome": "17",
                "mapped_hgvs_g_start": "100",
                "mapped_hgvs_g_stop": "100",
                "mapped_hgvs_g_ref": "A",
                "mapped_hgvs_g_alt": "G",
                "mapped_hgvs_p": "NP_999999.1:p.Cys4Arg",
            },
            {
                "gene_symbol": "TSC2",
                "mapped_hgvs_g": "NC_000016.10:g.999A>T",
                "mapped_hgvs_g_chromosome": "16",
                "mapped_hgvs_g_start": "999",
                "mapped_hgvs_g_stop": "999",
                "mapped_hgvs_g_ref": "A",
                "mapped_hgvs_g_alt": "T",
                "mapped_hgvs_p": "NP_001071.3:p.Asn1643His",
            },
        ],
        [
            "gene_symbol", "mapped_hgvs_g", "mapped_hgvs_g_chromosome", "mapped_hgvs_g_start",
            "mapped_hgvs_g_stop", "mapped_hgvs_g_ref", "mapped_hgvs_g_alt", "mapped_hgvs_p",
        ],
    )

    mod.main([
        str(in_path), str(out_path),
        "--revel-training-file", str(revel_training_path),
        "--mutpred2-training-file", str(mutpred2_training_path),
    ])

    rows = _read_tsv(out_path)
    assert rows[0]["revel.train"] == "true"
    assert rows[0]["mutpred2.train"] == "false"
    assert rows[1]["revel.train"] == "false"
    assert rows[1]["mutpred2.train"] == "true"


def test_main_requires_at_least_one_source_including_training_files(tmp_path):
    in_path = tmp_path / "in.tsv"
    out_path = tmp_path / "out.tsv"
    _write_tsv(in_path, [{"mapped_hgvs_g": "NC_000001.11:g.69094T>A"}], ["mapped_hgvs_g"])

    with pytest.raises(SystemExit):
        mod.main([str(in_path), str(out_path)])
