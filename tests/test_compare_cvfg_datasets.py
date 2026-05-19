import pytest
import pandas as pd

import src.compare_cvfg_datasets as mod


# ---------------------------------------------------------------------------
# _hgvs_c_expected_alt
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("hgvs_c, expected_alt", [
    # substitutions
    ("c.100A>G",          "G"),
    ("c.100_102AGC>TGA",  "TGA"),
    ("n.50C>T",           "T"),
    # delins
    ("c.235_237delinsCGC", "CGC"),
    ("c.100delinsTT",      "TT"),
    # deletions -> empty string
    ("c.100del",           ""),
    ("c.100_102del",       ""),
    ("c.100delA",          ""),
    # insertions
    ("c.100_101insAAT",    "AAT"),
    # duplications with explicit sequence
    ("c.100dupA",          "AA"),
    ("c.100_102dupAGC",    "AGCAGC"),
    # unresolvable cases -> None
    ("c.100dup",           None),   # no explicit sequence
    ("c.100_105inv",       None),   # inversion
    ("p.Gly12Val",         None),   # protein HGVS
    ("",                   None),
])
def test_hgvs_c_expected_alt(hgvs_c, expected_alt):
    assert mod._hgvs_c_expected_alt(hgvs_c) == expected_alt


# ---------------------------------------------------------------------------
# _is_transcript_alt_mismatch
# ---------------------------------------------------------------------------

def _row(hgvs_c, transcript_alt):
    return pd.Series({"hgvs_c": hgvs_c, "transcript_alt": transcript_alt})


@pytest.mark.parametrize("hgvs_c, transcript_alt, expect_error", [
    # the motivating example
    ("c.235_237delinsCGC",  "C",    True),
    ("c.235_237delinsCGC",  "CGC",  False),
    # substitutions
    ("c.100A>G",            "G",    False),
    ("c.100A>G",            "T",    True),
    ("c.100_102AGC>TGA",    "TGA",  False),
    ("c.100_102AGC>TGA",    "TG",   True),
    # case-insensitive
    ("c.100A>G",            "g",    False),
    # deletions
    ("c.100_102del",        "",     False),
    ("c.100_102del",        "A",    True),
    # insertions
    ("c.100_101insAAT",     "AAT",  False),
    ("c.100_101insAAT",     "AA",   True),
    # duplications with explicit sequence
    ("c.100dupA",           "AA",   False),
    # unresolvable -> always False (no false positives)
    ("c.100dup",            "X",    False),
    ("c.100_105inv",        "X",    False),
    # blank fields -> always False
    ("",                    "CGC",  False),
    ("c.235_237delinsCGC",  "",     False),
])
def test_is_transcript_alt_mismatch(hgvs_c, transcript_alt, expect_error):
    assert mod._is_transcript_alt_mismatch(_row(hgvs_c, transcript_alt)) == expect_error


# ---------------------------------------------------------------------------
# _is_vcf_hgvs_c_length_mismatch  (pre-join error filter)
# ---------------------------------------------------------------------------

def _vcf_row(ref, alt, hgvs_c):
    return pd.Series({"ref_allele": ref, "alt_allele": alt, "hgvs_c": hgvs_c})


@pytest.mark.parametrize("ref, alt, hgvs_c, expect_error", [
    # the motivating example: 2-bp VCF but 1-bp HGVS
    ("CT", "GG", "c.5570A>C",         True),
    # correct SNV
    ("A",  "C",  "c.5570A>C",         False),
    # correct 2-base substitution
    ("AG", "CC", "c.5569_5570AG>CC",  False),
    # 2-bp VCF alleles with matching 2-bp HGVS (minus strand)
    ("CT", "GG", "c.5569_5570AG>CC",  False),
    # non-substitution HGVS -> skip
    ("A",  "C",  "c.5570del",         False),
    ("TC", "T",  "c.5570delA",        False),
    # blank fields -> skip
    ("",   "C",  "c.5570A>C",         False),
    ("A",  "C",  "",                  False),
])
def test_is_vcf_hgvs_c_length_mismatch(ref, alt, hgvs_c, expect_error):
    assert mod._is_vcf_hgvs_c_length_mismatch(_vcf_row(ref, alt, hgvs_c)) == expect_error


# ---------------------------------------------------------------------------
# _merge_relaxed_start_off
# ---------------------------------------------------------------------------

def _make_a(rows):
    return pd.DataFrame(rows)


def _make_b(rows):
    return pd.DataFrame(rows)


def test_merge_relaxed_start_off_matches_off_by_one():
    df_a = _make_a([
        # correct span (3) == len(ref) (3) -> should NOT qualify
        {"Dataset": "DS1", "hg38_start": 10, "hg38_end": 12, "ref_allele": "AGC",
         "alt_allele": "TGA", "strand": 1, "__row_idx_a__": 0},
        # span=3 != len(ref)=2, but end-(start+1)+1=2 == len(ref) -> off by one, plus strand
        {"Dataset": "DS1", "hg38_start": 10, "hg38_end": 12, "ref_allele": "AG",
         "alt_allele": "TC", "strand": 1, "__row_idx_a__": 1},
        # off by one, minus strand: alt should match revcomp(c_alt)
        {"Dataset": "DS1", "hg38_start": 20, "hg38_end": 22, "ref_allele": "AG",
         "alt_allele": "CT", "strand": -1, "__row_idx_a__": 2},
    ])
    df_b = _make_b([
        # matches row 1 (adj_start=11): ref=AG, c_alt=TC
        {"dataset_name": "DS1", "mapped_hgvs_g_start": 11, "mapped_hgvs_g_stop": 12,
         "mapped_hgvs_g_ref": "AG", "mapped_hgvs_g_alt": "XX",
         "mapped_hgvs_c_alt": "TC", "__row_idx_b__": 0},
        # matches row 2 (adj_start=21): ref=AG, revcomp(c_alt=AG)=CT
        {"dataset_name": "DS1", "mapped_hgvs_g_start": 21, "mapped_hgvs_g_stop": 22,
         "mapped_hgvs_g_ref": "AG", "mapped_hgvs_g_alt": "YY",
         "mapped_hgvs_c_alt": "AG", "__row_idx_b__": 1},
    ])
    result = mod._merge_relaxed_start_off(df_a, df_b)
    assert sorted(result["__row_idx_a__"].tolist()) == [1, 2]
    assert sorted(result["__row_idx_b__"].tolist()) == [0, 1]


def test_merge_relaxed_start_off_no_match_when_span_correct():
    """Rows whose span already matches len(ref) are not touched."""
    df_a = _make_a([
        {"Dataset": "DS1", "hg38_start": 10, "hg38_end": 12, "ref_allele": "AGC",
         "alt_allele": "TGA", "strand": 1, "__row_idx_a__": 0},
    ])
    df_b = _make_b([
        {"dataset_name": "DS1", "mapped_hgvs_g_start": 11, "mapped_hgvs_g_stop": 12,
         "mapped_hgvs_g_ref": "AGC", "mapped_hgvs_g_alt": "XX",
         "mapped_hgvs_c_alt": "TGA", "__row_idx_b__": 0},
    ])
    result = mod._merge_relaxed_start_off(df_a, df_b)
    assert result.empty


def test_merge_relaxed_start_off_wrong_transcript_alt():
    """Off-by-one start row is not matched when transcript alt doesn't agree."""
    df_a = _make_a([
        {"Dataset": "DS1", "hg38_start": 10, "hg38_end": 12, "ref_allele": "AG",
         "alt_allele": "TC", "strand": 1, "__row_idx_a__": 0},
    ])
    df_b = _make_b([
        {"dataset_name": "DS1", "mapped_hgvs_g_start": 11, "mapped_hgvs_g_stop": 12,
         "mapped_hgvs_g_ref": "AG", "mapped_hgvs_g_alt": "XX",
         "mapped_hgvs_c_alt": "GG",  # wrong
         "__row_idx_b__": 0},
    ])
    result = mod._merge_relaxed_start_off(df_a, df_b)
    assert result.empty
