#!/usr/bin/env python3
"""
Compare integrated_variant_effect_dataset_condensed.tsv (A) with
a CVFG variants flat file (B) via a positional-genome join.

Join key (compound, strict):
  A.Dataset          == B.dataset_name
  A.hg38_start       == B.mapped_hgvs_g_start
  A.hg38_end         == B.mapped_hgvs_g_stop
  A.ref_allele       == B.mapped_hgvs_g_ref
  A.alt_allele       == B.mapped_hgvs_g_alt

Relaxed fallback join 1 (for intronic codon-split variants):
  Same as strict join except the alt check uses B.mapped_hgvs_c_alt (the
  transcript-level alt, intron-free) instead of B.mapped_hgvs_g_alt.
  On the minus strand (strand == -1) A's genomic alt is the reverse complement
  of the transcript alt.

Relaxed fallback join 2 (for VCF-anchored deletions):
  File A stores deletions with an extra anchor nucleotide included in ref and
  alt, following VCF convention (e.g. ref=TC alt=C or ref=AA alt=A).
  File B stores only the deleted bases with a blank alt.
  Three position conventions are handled for A:
    Case 1 (span == len(ref)):  anchor included in span, deletion starts after it.
    Case 2 (span == del_len):   span covers only the deleted bases.
    Case 3 (stop blank):        only start is given (= first deleted position).
  Both left-anchor (ref.startswith(alt)) and right-anchor (ref.endswith(alt))
  orientations are tried.

Relaxed fallback join 3 (for off-by-one start coordinate in A):
  Some A rows have len(ref) != end - start + 1; trying start+1 and matching
  against B's transcript alt (with strand-aware revcomp) resolves them.

Pre-filter (applied before the join):
  compare_errors_in_a.tsv         - A rows with VCF/HGVS-c allele-length mismatch

Outputs (written to --output-dir):
  compare_only_in_a.tsv           - clean A rows with no match in B (all A columns)
  compare_only_in_a_coord_errors.tsv - unmatched A rows with coordinate or ref-allele errors
  compare_only_in_b.tsv           - rows from B with no match in A (all B columns)
  compare_joined.tsv              - matched rows with selected columns from both
"""

import argparse
import os
import re
import sys
from functools import lru_cache
from typing import Optional, Any

import pandas as pd

# Matches c./n. substitutions: c.123A>C or c.123_125AGC>TTA
_HGVS_C_SUBST_RE = re.compile(
    r"^[cn]\.[0-9*+\-]+(?:_[0-9*+\-]+)?([A-Za-z]+)>([A-Za-z]+)$"
)
# Matches c./n. delins: c.235_237delinsCGC or c.235delinsCGC
_HGVS_C_DELINS_RE = re.compile(
    r"^[cn]\.[0-9*+\-]+(?:_[0-9*+\-]+)?delins([A-Za-z]+)$"
)
# Matches c./n. insertions: c.235_236insCGC
_HGVS_C_INS_RE = re.compile(
    r"^[cn]\.[0-9*+\-]+_[0-9*+\-]+ins([A-Za-z]+)$"
)
# Matches c./n. duplications: c.235dup or c.235_237dup
_HGVS_C_DUP_RE = re.compile(
    r"^[cn]\.[0-9*+\-]+(?:_[0-9*+\-]+)?dup([A-Za-z]*)$"
)


DEFAULT_FILE_A = "data/cvfg/integrated_variant_effect_dataset.tsv"
DEFAULT_FILE_B = "data/cvfg/v5/cvfg_variants.3.flat.tsv"
DEFAULT_FILE_C = "data/cvfg/dataset_to_score_set.tsv"
DEFAULT_OUTPUT_DIR = "data/cvfg"

SPLICEAI_COLS = ["spliceAI_DS_AG", "spliceAI_DS_AL", "spliceAI_DS_DG", "spliceAI_DS_DL"]

# Ordered columns for the joined output TSV
JOINED_OUTPUT_COLS = [
    "ID",
    "dataset_name",
    "gene_symbol",
    "variant_urn",
    "raw_hgvs_nt",
    "raw_hgvs_pro",
    "mapped_hgvs_g",
    "mapped_hgvs_c",
    "mapped_hgvs_p",
    "mapping_error",
    "reverse_translation_error",
    "clingen_allele_id",
    "dna_clingen_allele_id",
    "score",
    "clinvar_sig_2018",
    "clinvar_sig_2025",
    "gnomad_MAF",
    "has_spliceai",
    "clinvar.201801.variation_id",
    "clinvar.201801.allele_id",
    "clinvar.201801.clinical_significance",
    "clinvar.202501.variation_id",
    "clinvar.202501.allele_id",
    "clinvar.202501.clinical_significance",
    "clinvar.202601.variation_id",
    "clinvar.202601.allele_id",
    "clinvar.202601.clinical_significance",
    "gnomad.v4_1.minor_allele_frequency",
    "spliceai.max_delta_score",
]

# GRCh38 chromosome label → RefSeq versioned accession (for UTA/seqrepo lookup)
_CHROMOSOME_TO_ACCESSION_HG38: dict[str, str] = {
    "1":  "NC_000001.11", "2":  "NC_000002.12", "3":  "NC_000003.12",
    "4":  "NC_000004.12", "5":  "NC_000005.10", "6":  "NC_000006.12",
    "7":  "NC_000007.14", "8":  "NC_000008.11", "9":  "NC_000009.12",
    "10": "NC_000010.11", "11": "NC_000011.10", "12": "NC_000012.12",
    "13": "NC_000013.11", "14": "NC_000014.9",  "15": "NC_000015.10",
    "16": "NC_000016.10", "17": "NC_000017.11", "18": "NC_000018.10",
    "19": "NC_000019.10", "20": "NC_000020.11", "21": "NC_000021.9",
    "22": "NC_000022.11", "X":  "NC_000023.11", "Y":  "NC_000024.10",
    "M":  "NC_012920.1",
}

# Module-level UTA data-provider handle (lazily initialised)
_UTA_HDP: Optional[Any] = None
_UTA_HDP_TRIED: bool = False


def _get_uta_hdp() -> Optional[Any]:
    """Return a connected UTA data provider, or None if UTA is unavailable."""
    global _UTA_HDP, _UTA_HDP_TRIED
    if _UTA_HDP_TRIED:
        return _UTA_HDP
    _UTA_HDP_TRIED = True
    uta_db_url = (os.environ.get("UTA_DB_URL") or "").strip()
    if not uta_db_url:
        return None
    try:
        import hgvs.dataproviders.uta  # type: ignore[import-untyped]
        _UTA_HDP = hgvs.dataproviders.uta.connect(uta_db_url)
    except Exception as exc:
        print(f"  Warning: UTA connection failed ({exc}); skipping genomic ref allele check")
    return _UTA_HDP


@lru_cache(maxsize=10000)
def _fetch_genomic_seq(accession: str, start_1based: int, stop_1based: int) -> Optional[str]:
    """Fetch the reference sequence for a 1-based inclusive genomic range from UTA."""
    hdp = _get_uta_hdp()
    if hdp is None:
        return None
    try:
        # UTA get_seq uses 0-based half-open intervals (Python slice convention)
        seq = hdp.get_seq(accession, start_1based - 1, stop_1based)
        return seq if seq else None
    except Exception:
        return None


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "file_a",
        nargs="?",
        default=DEFAULT_FILE_A,
        help="Integrated variant effect dataset (condensed TSV). Default: %(default)s",
    )
    parser.add_argument(
        "file_b",
        nargs="?",
        default=DEFAULT_FILE_B,
        help="CVFG variants with scores (flat TSV). Default: %(default)s",
    )
    parser.add_argument(
        "file_c",
        nargs="?",
        default=DEFAULT_FILE_C,
        help="Dataset-to-score-set URN mapping TSV. Default: %(default)s",
    )
    parser.add_argument(
        "--output-dir",
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for output TSV files. Default: %(default)s",
    )
    parser.add_argument(
        "--check-ref-alleles",
        action=argparse.BooleanOptionalAction,
        default=True,
        help=(
            "Verify that ref_allele matches the actual GRCh38 sequence at the specified "
            "hg38 coordinates (requires UTA_DB_URL). Mismatching rows are moved to "
            "compare_only_in_a_coord_errors.tsv. Default: %(default)s"
        ),
    )
    return parser.parse_args()


def read_tsv(path: str, label: str) -> pd.DataFrame:
    print(f"Reading {label}: {path}")
    if not os.path.isfile(path):
        sys.exit(f"Error: file not found: {path}")
    return pd.read_csv(path, sep="\t", low_memory=False)


def _hgvs_c_expected_alt(hgvs_c: str) -> str | None:
    """Return the expected transcript alt nucleotides for a c./n. HGVS string.

    Returns:
        A DNA string (upper-case) when the alt can be inferred unambiguously:
          - substitution / delins  → the inserted/alt sequence
          - insertion              → the inserted sequence
          - duplication            → the duplicated sequence (when explicit)
          - deletion               → empty string ""
        Returns None when the operation type is unrecognised or unresolvable
        (e.g. a dup without an explicit sequence).
    """
    s = (hgvs_c or "").strip()
    if not s:
        return None
    m = _HGVS_C_SUBST_RE.match(s)
    if m:
        return m.group(2).upper()
    m = _HGVS_C_DELINS_RE.match(s)
    if m:
        return m.group(1).upper()
    m = _HGVS_C_INS_RE.match(s)
    if m:
        return m.group(1).upper()
    m = _HGVS_C_DUP_RE.match(s)
    if m:
        dup_seq = m.group(1)
        return (dup_seq * 2).upper() if dup_seq else None
    if re.match(r"^[cn]\.[0-9*+\-]+(?:_[0-9*+\-]+)?del[A-Za-z]*$", s):
        return ""
    return None


def _is_transcript_alt_mismatch(row: pd.Series) -> bool:
    """True when transcript_alt disagrees with what hgvs_c implies.

    Only fires when both columns are non-blank and the expected alt can be
    inferred.  Comparison is case-insensitive.
    """
    hgvs_c = str(row.get("hgvs_c") or "").strip()
    transcript_alt = str(row.get("transcript_alt") or "").strip()
    if not hgvs_c or not transcript_alt:
        return False
    expected = _hgvs_c_expected_alt(hgvs_c)
    if expected is None:
        return False
    return transcript_alt.upper() != expected


def _hgvs_c_subst_lengths(hgvs_c: str) -> tuple[int, int] | None:
    """Return (len_ref, len_alt) for a c./n. substitution HGVS string, or None."""
    m = _HGVS_C_SUBST_RE.match((hgvs_c or "").strip())
    if not m:
        return None
    return len(m.group(1)), len(m.group(2))


def _is_vcf_hgvs_c_length_mismatch(row: pd.Series) -> bool:
    """Return True if VCF allele lengths disagree with the HGVS c. substitution.

    For substitution variants (hgvs_c contains '>'), the number of nucleotides
    in ref_allele / alt_allele must equal the number in the HGVS c. ref / alt.
    Length is strand-invariant, so no strand flip is needed.
    """
    hgvs_c = str(row.get("hgvs_c") or "").strip()
    if not hgvs_c:
        return False
    lengths = _hgvs_c_subst_lengths(hgvs_c)
    if lengths is None:
        return False  # not a substitution, or couldn't parse
    c_ref_len, c_alt_len = lengths
    ref = str(row.get("ref_allele") or "").strip()
    alt = str(row.get("alt_allele") or "").strip()
    if not ref or not alt:
        return False
    return len(ref) != c_ref_len or len(alt) != c_alt_len


def _is_genomic_ref_allele_mismatch(row: pd.Series) -> bool:
    """True when ref_allele doesn't match the actual GRCh38 sequence at hg38 coordinates.

    Uses UTA to look up the reference sequence.  Returns False (no error flagged) when
    the UTA resolver is unavailable, or when any required field is missing/unparseable.
    Comparison is case-insensitive.
    """
    chrom = str(row.get("Chrom") or "").strip().lstrip("chr")
    ref = str(row.get("ref_allele") or "").strip()
    if not chrom or not ref:
        return False
    accession = _CHROMOSOME_TO_ACCESSION_HG38.get(chrom)
    if not accession:
        return False
    try:
        start = int(row["hg38_start"])
        stop = int(row["hg38_end"])
    except (TypeError, ValueError):
        return False
    actual = _fetch_genomic_seq(accession, start, stop)
    if actual is None:
        return False  # resolver unavailable; don't falsely flag
    return actual.upper() != ref.upper()


def add_has_spliceai(df: pd.DataFrame) -> pd.DataFrame:
    """Add has_spliceai: True if any spliceAI DS column is non-null and non-empty."""
    present = [c for c in SPLICEAI_COLS if c in df.columns]
    df = df.copy()
    if not present:
        df["has_spliceai"] = False
        return df
    def _any_nonempty(row: pd.Series) -> bool:
        return any(pd.notna(row[c]) and str(row[c]).strip() != "" for c in present)
    df["has_spliceai"] = df.apply(_any_nonempty, axis=1)
    return df


def coerce_int_key(series: pd.Series) -> pd.Series:
    """Convert a genomic position column to nullable Int64 for safe joining."""
    return pd.to_numeric(series, errors="coerce").astype("Int64")


def recover_original_cols(
    merged: pd.DataFrame,
    original_cols: list[str],
    suffix: str,
) -> pd.DataFrame:
    """
    From a merged DataFrame, reconstruct a DataFrame with the original column names.
    Prefers the suffixed version (e.g. 'col_a') over the plain name when both exist.
    """
    col_mapping: list[tuple[str, str]] = []  # (merged_col, original_col)
    for col in original_cols:
        suffixed = col + suffix
        if suffixed in merged.columns:
            col_mapping.append((suffixed, col))
        elif col in merged.columns:
            col_mapping.append((col, col))
        # If neither exists (shouldn't happen), skip
    merged_names = [mc for mc, _ in col_mapping]
    original_names = [orig for _, orig in col_mapping]
    result = merged[merged_names].copy()
    result.columns = original_names
    return result


def resolve_col(merged: pd.DataFrame, name: str, prefer_suffix: str) -> str | None:
    """Return the actual column name in `merged` for the logical column `name`."""
    suffixed = name + prefer_suffix
    if suffixed in merged.columns:
        return suffixed
    if name in merged.columns:
        return name
    return None


def _revcomp(seq: str) -> str:
    """Return the reverse complement of a DNA sequence string."""
    comp = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp)[::-1]


def _merge_relaxed_start_off(sub_a: pd.DataFrame, sub_b: pd.DataFrame) -> pd.DataFrame:
    """
    Inner merge for rows where A's hg38_start is off by +1 relative to the
    ref allele length (i.e. len(ref) != end - start + 1 but
    len(ref) == end - (start+1) + 1 == end - start).

    We rebuild the join key using hg38_start + 1 and match alt using
    B's transcript alt (mapped_hgvs_c_alt) with strand-aware revcomp,
    identical to _merge_relaxed_alt.
    """
    # Only rows where start is off by 1 in the expected direction.
    def _is_start_off(row) -> bool:
        try:
            start = int(row["hg38_start"])
            stop = int(row["hg38_end"])
        except (TypeError, ValueError):
            return False
        ref = str(row.get("ref_allele") or "")
        if not ref:
            return False
        span = stop - start + 1
        return span != len(ref) and (stop - (start + 1) + 1) == len(ref)

    cand_a = sub_a[sub_a.apply(_is_start_off, axis=1)].copy()
    if cand_a.empty:
        return pd.DataFrame()

    cand_a = cand_a.copy()
    cand_a = cand_a.assign(_adj_start=cand_a["hg38_start"] + 1)

    candidate = cand_a.merge(
        sub_b,
        left_on=["Dataset", "_adj_start", "hg38_end", "ref_allele"],
        right_on=[
            "dataset_name",
            "mapped_hgvs_g_start",
            "mapped_hgvs_g_stop",
            "mapped_hgvs_g_ref",
        ],
        how="inner",
        suffixes=("_a", "_b"),
    )
    if candidate.empty:
        return pd.DataFrame()

    def _matches(row) -> bool:
        a_alt = row.get("alt_allele")
        c_alt = row.get("mapped_hgvs_c_alt")
        if pd.isna(a_alt) or pd.isna(c_alt):
            return False
        a_alt, c_alt = str(a_alt), str(c_alt)
        if not a_alt or not c_alt:
            return False
        try:
            strand_val = int(row.get("strand"))
        except (TypeError, ValueError):
            strand_val = 1
        return a_alt == (_revcomp(c_alt) if strand_val == -1 else c_alt)

    result = candidate[candidate.apply(_matches, axis=1)].copy()
    return result.drop(columns=["_adj_start"], errors="ignore")


def _merge_relaxed_alt(sub_a: pd.DataFrame, sub_b: pd.DataFrame) -> pd.DataFrame:
    """
    Inner merge on position + ref only (no alt), then keep rows where A's
    alt_allele equals B's transcript alt (mapped_hgvs_c_alt), accounting for
    strand orientation.  On the minus strand (strand == -1) the genomic alt
    allele stored in A is the reverse complement of the transcript alt, so
    the comparison is A.alt_allele == revcomp(B.mapped_hgvs_c_alt).
    """
    candidate = sub_a.merge(
        sub_b,
        left_on=["Dataset", "hg38_start", "hg38_end", "ref_allele"],
        right_on=[
            "dataset_name",
            "mapped_hgvs_g_start",
            "mapped_hgvs_g_stop",
            "mapped_hgvs_g_ref",
        ],
        how="inner",
        suffixes=("_a", "_b"),
    )
    if candidate.empty:
        return candidate

    def _matches(row) -> bool:
        a_alt = row.get("alt_allele")
        c_alt = row.get("mapped_hgvs_c_alt")
        if pd.isna(a_alt) or pd.isna(c_alt):
            return False
        a_alt, c_alt = str(a_alt), str(c_alt)
        if not a_alt or not c_alt:
            return False
        try:
            strand_val = int(row.get("strand"))
        except (TypeError, ValueError):
            strand_val = 1  # assume plus strand when unknown
        return a_alt == (_revcomp(c_alt) if strand_val == -1 else c_alt)

    return candidate[candidate.apply(_matches, axis=1)].copy()


def _merge_relaxed_del(sub_a: pd.DataFrame, sub_b: pd.DataFrame) -> pd.DataFrame:
    """
    Inner merge for deletion variants where A uses VCF-style anchored format
    (one extra anchor nucleotide included in ref/alt) but B stores only the
    deleted bases with a blank alt.

    For each A row detected as a VCF-anchored deletion (ref longer than alt,
    and alt is a prefix or suffix of ref), candidate (b_start, b_stop,
    del_bases) tuples are computed for all three position conventions and both
    anchor orientations, then matched against B rows that have a blank alt.
    """

    def _del_candidates(row) -> list[tuple[int, int, str]]:
        ref = str(row.get("ref_allele") or "")
        alt = str(row.get("alt_allele") or "")
        if not alt or len(ref) <= len(alt):
            return []
        try:
            a_start = int(row["hg38_start"])
        except (TypeError, ValueError):
            return []
        a_stop_raw = row["hg38_end"]
        a_stop_blank = pd.isna(a_stop_raw)
        try:
            a_stop = None if a_stop_blank else int(a_stop_raw)
        except (TypeError, ValueError):
            return []

        cands: set[tuple[int, int, str]] = set()

        def _add_cands(del_bases: str, is_left: bool) -> None:
            if not del_bases:
                return
            del_len = len(del_bases)
            if a_stop_blank:  # case 3: only start given
                cands.add((a_start, a_start + del_len - 1, del_bases))
            else:
                span = a_stop - a_start + 1  # type: ignore[operator]
                if span == len(ref):  # case 1: span includes anchor
                    if is_left:  # anchor at left → deletion starts after anchor
                        cands.add((a_start + len(alt), a_stop, del_bases))
                    else:  # anchor at right → deletion ends before anchor
                        cands.add((a_start, a_stop - len(alt), del_bases))
                if span == del_len:  # case 2: span covers deletion only
                    cands.add((a_start, a_stop, del_bases))

        if ref.startswith(alt):
            _add_cands(ref[len(alt):], is_left=True)
        if ref.endswith(alt) and not ref.startswith(alt):
            _add_cands(ref[: -len(alt)], is_left=False)
        # Homopolymer (startswith AND endswith): case 1 differs by anchor side;
        # handle both positions explicitly.
        if ref.startswith(alt) and ref.endswith(alt) and not a_stop_blank:
            span = a_stop - a_start + 1  # type: ignore[operator]
            if span == len(ref):
                del_l = ref[len(alt):]
                del_r = ref[: -len(alt)]
                if del_l != del_r or a_start + len(alt) != a_start:
                    _add_cands(del_r, is_left=False)  # right-anchor case 1

        return list(cands)

    # Build exploded A-side: one row per candidate key.
    a_exp_rows: list[dict] = []
    for _, row in sub_a.iterrows():
        for b_start, b_stop, del_bases in _del_candidates(row):
            a_exp_rows.append({
                **row.to_dict(),
                "_dbk": "%d:%d:%s" % (b_start, b_stop, del_bases),
            })
    if not a_exp_rows:
        return pd.DataFrame()
    a_exp = pd.DataFrame(a_exp_rows)

    # Filter B to deletion rows (alt blank, ref non-empty).
    alt_blank = sub_b["mapped_hgvs_g_alt"].isna() | (
        sub_b["mapped_hgvs_g_alt"].astype(str).str.strip().isin(["", "nan"])
    )
    ref_nonempty = sub_b["mapped_hgvs_g_ref"].notna() & ~(
        sub_b["mapped_hgvs_g_ref"].astype(str).str.strip().isin(["", "nan"])
    )
    sub_b_del = sub_b[alt_blank & ref_nonempty].copy()
    if sub_b_del.empty:
        return pd.DataFrame()

    sub_b_del["_dbk"] = (
        sub_b_del["mapped_hgvs_g_start"].astype(str)
        + ":"
        + sub_b_del["mapped_hgvs_g_stop"].astype(str)
        + ":"
        + sub_b_del["mapped_hgvs_g_ref"].astype(str)
    )

    result = a_exp.merge(
        sub_b_del,
        left_on=["Dataset", "_dbk"],
        right_on=["dataset_name", "_dbk"],
        how="inner",
        suffixes=("_a", "_b"),
    )
    if result.empty:
        return result

    result = result.drop_duplicates(subset=["__row_idx_a__", "__row_idx_b__"])
    return result.drop(columns=["_dbk", "_dbk_a", "_dbk_b"], errors="ignore")


def main() -> None:
    args = parse_args()

    # ── 1. Read inputs ──────────────────────────────────────────────────────────
    df_a = read_tsv(args.file_a, "A (condensed dataset)")
    df_b = read_tsv(args.file_b, "B (cvfg variants)")
    df_c = read_tsv(args.file_c, "C (dataset-to-URN mapping)")
    print(f"  A: {len(df_a):,} rows, {len(df_a.columns)} columns")
    print(f"  B: {len(df_b):,} rows, {len(df_b.columns)} columns")
    print(f"  C: {len(df_c):,} rows")

    # ── 2. Add score_set_urn to B ───────────────────────────────────────────────
    df_b = df_b.copy()
    df_b["score_set_urn"] = df_b["variant_urn"].str.split("#").str[0]

    # ── 3. Add dataset_name to B via C ─────────────────────────────────────────
    urn_to_dataset: dict[str, str] = df_c.set_index("score_set_urn")["dataset_name"].to_dict()
    df_b["dataset_name"] = df_b["score_set_urn"].map(urn_to_dataset)
    unmapped = df_b["dataset_name"].isna().sum()
    if unmapped:
        print(f"  Warning: {unmapped:,} rows in B have no dataset_name mapping in C")

    # ── 4. Add has_spliceai to A ────────────────────────────────────────────────
    df_a = add_has_spliceai(df_a)

    # ── 4b. Filter erroneous rows from A (VCF/HGVS-c length mismatch) ──────────
    error_mask = df_a.apply(_is_vcf_hgvs_c_length_mismatch, axis=1)
    df_a_errors = df_a[error_mask].copy()
    df_a = df_a[~error_mask].copy()
    out_errors_a = os.path.join(args.output_dir, "compare_errors_in_a.tsv")
    os.makedirs(args.output_dir, exist_ok=True)
    df_a_errors.to_csv(out_errors_a, sep="\t", index=False)
    print(f"  Erroneous rows in A (VCF/HGVS-c mismatch): {len(df_a_errors):,} → {out_errors_a}")
    print(f"  Clean rows in A for join: {len(df_a):,}")

    # ── 5. Normalize join-key types ─────────────────────────────────────────────
    # A stores positions as floats (e.g. 32356440.0); B stores them as integers.
    df_a = df_a.copy()
    df_a["hg38_start"] = coerce_int_key(df_a["hg38_start"])
    df_a["hg38_end"] = coerce_int_key(df_a["hg38_end"])
    df_b["mapped_hgvs_g_start"] = coerce_int_key(df_b["mapped_hgvs_g_start"])
    df_b["mapped_hgvs_g_stop"] = coerce_int_key(df_b["mapped_hgvs_g_stop"])

    # Save original column lists (before tracking indices are added) so that
    # recover_original_cols output is not polluted by the tracking columns.
    a_orig_cols = list(df_a.columns)
    b_orig_cols = list(df_b.columns)

    # Add row-tracking indices to support the relaxed join below.
    df_a["__row_idx_a__"] = range(len(df_a))
    df_b["__row_idx_b__"] = range(len(df_b))

    # ── 6. Outer merge with indicator ───────────────────────────────────────────
    merged = df_a.merge(
        df_b,
        left_on=["Dataset", "hg38_start", "hg38_end", "ref_allele", "alt_allele"],
        right_on=[
            "dataset_name",
            "mapped_hgvs_g_start",
            "mapped_hgvs_g_stop",
            "mapped_hgvs_g_ref",
            "mapped_hgvs_g_alt",
        ],
        how="outer",
        indicator=True,
        suffixes=("_a", "_b"),
    )

    only_a_mask = merged["_merge"] == "left_only"
    only_b_mask = merged["_merge"] == "right_only"
    both_mask = merged["_merge"] == "both"

    print(f"\nJoin results (strict):")
    print(f"  Rows only in A:  {only_a_mask.sum():,}")
    print(f"  Rows only in B:  {only_b_mask.sum():,}")
    print(f"  Rows matched:    {both_mask.sum():,}")

    # ── 6b. Relaxed join for rows unmatched due to alt-allele representation ────
    # For codon changes spanning an intron, file A stores only the substituted
    # codon as alt_allele while file B stores the full genomic alt (codon + intron).
    # We match these by comparing A.alt_allele to B.mapped_hgvs_c_alt (the
    # transcript-level alt, intron-free).  For minus-strand genes (strand == -1)
    # A's alt_allele is the reverse complement of B's transcript alt.
    unmatched_a_idxs = set(
        merged.loc[only_a_mask, "__row_idx_a__"].dropna().astype(int)
    )
    unmatched_b_idxs = set(
        merged.loc[only_b_mask, "__row_idx_b__"].dropna().astype(int)
    )
    sub_a = df_a[df_a["__row_idx_a__"].isin(unmatched_a_idxs)].copy()
    sub_b = df_b[df_b["__row_idx_b__"].isin(unmatched_b_idxs)].copy()
    relaxed_both = _merge_relaxed_alt(sub_a, sub_b)

    relaxed_matched_a_idxs: set[int] = (
        set(relaxed_both["__row_idx_a__"].astype(int))
        if not relaxed_both.empty
        else set()
    )
    relaxed_matched_b_idxs: set[int] = (
        set(relaxed_both["__row_idx_b__"].astype(int))
        if not relaxed_both.empty
        else set()
    )

    # ── 6c. Second relaxed join: VCF-anchored deletion discrepancy ──────────────
    # A stores deletions with an anchor nucleotide in ref/alt (VCF convention);
    # B stores only the deleted bases with a blank alt.  Run on rows that are
    # still unmatched after the strict and first relaxed joins.
    sub_a2 = df_a[
        df_a["__row_idx_a__"].isin(unmatched_a_idxs - relaxed_matched_a_idxs)
    ].copy()
    sub_b2 = df_b[
        df_b["__row_idx_b__"].isin(unmatched_b_idxs - relaxed_matched_b_idxs)
    ].copy()
    relaxed_del_both = _merge_relaxed_del(sub_a2, sub_b2)

    relaxed_del_matched_a_idxs: set[int] = (
        set(relaxed_del_both["__row_idx_a__"].astype(int))
        if not relaxed_del_both.empty
        else set()
    )
    relaxed_del_matched_b_idxs: set[int] = (
        set(relaxed_del_both["__row_idx_b__"].astype(int))
        if not relaxed_del_both.empty
        else set()
    )

    # ── 6d. Third relaxed join: off-by-one start coordinate + transcript alt ────
    # Some A rows have len(ref) != end - start + 1; trying start+1 and matching
    # against B's transcript alt (with strand-aware revcomp) resolves them.
    sub_a3 = df_a[
        df_a["__row_idx_a__"].isin(
            unmatched_a_idxs - relaxed_matched_a_idxs - relaxed_del_matched_a_idxs
        )
    ].copy()
    sub_b3 = df_b[
        df_b["__row_idx_b__"].isin(
            unmatched_b_idxs - relaxed_matched_b_idxs - relaxed_del_matched_b_idxs
        )
    ].copy()
    relaxed_startoff_both = _merge_relaxed_start_off(sub_a3, sub_b3)

    relaxed_startoff_matched_a_idxs: set[int] = (
        set(relaxed_startoff_both["__row_idx_a__"].astype(int))
        if not relaxed_startoff_both.empty
        else set()
    )
    relaxed_startoff_matched_b_idxs: set[int] = (
        set(relaxed_startoff_both["__row_idx_b__"].astype(int))
        if not relaxed_startoff_both.empty
        else set()
    )

    all_relaxed_a = relaxed_matched_a_idxs | relaxed_del_matched_a_idxs | relaxed_startoff_matched_a_idxs
    all_relaxed_b = relaxed_matched_b_idxs | relaxed_del_matched_b_idxs | relaxed_startoff_matched_b_idxs

    truly_only_a_mask = only_a_mask & ~merged["__row_idx_a__"].isin(all_relaxed_a)
    truly_only_b_mask = only_b_mask & ~merged["__row_idx_b__"].isin(all_relaxed_b)

    if relaxed_matched_a_idxs:
        print(f"  Relaxed-matched (transcript-alt discrepancy): {len(relaxed_matched_a_idxs):,}")
    if relaxed_del_matched_a_idxs:
        print(f"  Relaxed-matched (deletion anchor discrepancy): {len(relaxed_del_matched_a_idxs):,}")
    if relaxed_startoff_matched_a_idxs:
        print(f"  Relaxed-matched (off-by-one start + transcript alt): {len(relaxed_startoff_matched_a_idxs):,}")
    if all_relaxed_a:
        print(f"  Adjusted rows only in A: {truly_only_a_mask.sum():,}")
        print(f"  Adjusted rows only in B: {truly_only_b_mask.sum():,}")

    # ── 7. Output: only-in-A (split into clean rows and coord/ref-length errors) ─
    only_a_df = recover_original_cols(
        merged.loc[truly_only_a_mask], a_orig_cols, "_a"
    )

    def _ref_coord_mismatch(row: pd.Series) -> bool:
        """True when len(ref_allele) != hg38_end - hg38_start + 1."""
        ref = str(row.get("ref_allele") or "").strip()
        if not ref:
            return False
        try:
            start = int(row["hg38_start"])
            stop = int(row["hg38_end"])
        except (TypeError, ValueError):
            return False
        return len(ref) != stop - start + 1

    coord_err_mask = only_a_df.apply(_ref_coord_mismatch, axis=1)
    transcript_alt_err_mask = only_a_df.apply(_is_transcript_alt_mismatch, axis=1)

    # Genomic ref-allele check (requires UTA; skipped when unavailable)
    if args.check_ref_alleles:
        hdp = _get_uta_hdp()
        if hdp is not None:
            print("\nChecking genomic ref alleles against GRCh38 sequence via UTA…")
            genomic_ref_err_mask = only_a_df.apply(_is_genomic_ref_allele_mismatch, axis=1)
            print(f"  Genomic ref-allele mismatches: {genomic_ref_err_mask.sum():,}")
        else:
            if args.check_ref_alleles:
                print(
                    "\n  Warning: UTA unavailable (set UTA_DB_URL); "
                    "skipping genomic ref allele check"
                )
            genomic_ref_err_mask = pd.Series(False, index=only_a_df.index)
    else:
        genomic_ref_err_mask = pd.Series(False, index=only_a_df.index)

    combined_err_mask = coord_err_mask | transcript_alt_err_mask | genomic_ref_err_mask
    only_a_coord_errors_df = only_a_df[combined_err_mask].copy()
    only_a_df = only_a_df[~combined_err_mask].copy()

    out_only_a_coord_errors = os.path.join(args.output_dir, "compare_only_in_a_coord_errors.tsv")
    only_a_coord_errors_df.to_csv(out_only_a_coord_errors, sep="\t", index=False)

    out_only_a = os.path.join(args.output_dir, "compare_only_in_a.tsv")
    only_a_df.to_csv(out_only_a, sep="\t", index=False)
    print(f"\nWrote only-in-A coord errors ({len(only_a_coord_errors_df):,} rows) → {out_only_a_coord_errors}")
    print(f"Wrote only-in-A ({len(only_a_df):,} rows) → {out_only_a}")

    # ── 8. Output: only-in-B ────────────────────────────────────────────────────
    only_b_df = recover_original_cols(
        merged.loc[truly_only_b_mask], b_orig_cols, "_b"
    )
    out_only_b = os.path.join(args.output_dir, "compare_only_in_b.tsv")
    only_b_df.to_csv(out_only_b, sep="\t", index=False)
    print(f"Wrote only-in-B ({len(only_b_df):,} rows) → {out_only_b}")

    # ── 9. Output: joined rows with selected columns ────────────────────────────

    # Columns sourced from A
    A_COLS = [
        "ID",
        "clinvar_sig_2018",
        "clinvar_sig_2025",
        "gnomad_MAF",
        "has_spliceai",
        "Uuid_ClinGen_repo",
        "auth_reported_func_class",
        "consequence",
        "AM_score",
        "MutPred2",
        "REVEL"
    ]
    # Columns sourced from B
    B_COLS = [
        "dataset_name",
        "gene_symbol",
        "variant_urn",
        "raw_hgvs_nt",
        "raw_hgvs_pro",
        "mapped_hgvs_g",
        "mapped_hgvs_c",
        "mapped_hgvs_p",
        "mapping_error",
        "reverse_translation_error",
        "clingen_allele_id",
        "dna_clingen_allele_id",
        "score",
        "clinvar.201801.variation_id",
        "clinvar.201801.allele_id",
        "clinvar.201801.clinical_significance",
        "clinvar.202501.variation_id",
        "clinvar.202501.allele_id",
        "clinvar.202501.clinical_significance",
        "clinvar.202601.variation_id",
        "clinvar.202601.allele_id",
        "clinvar.202601.clinical_significance",
        "gnomad.v4_1.minor_allele_frequency",
        "spliceai.max_delta_score",
        "clingen_evidence_repository.Uuid",
        "mavedb.primary_calibration.functional_class",
        "vep.mutational_consequence",
        "alphamissense.pathogenicity",
        "mutpred2.score",
        "revel.score"
    ]

    def _build_joined_df(rows_df: pd.DataFrame) -> pd.DataFrame:
        """Extract A and B columns from a merged DataFrame into the joined output shape."""
        out = pd.DataFrame()
        for logical in A_COLS:
            actual = resolve_col(rows_df, logical, "_a")
            out[logical] = rows_df[actual].values if actual else pd.NA
        for logical in B_COLS:
            actual = resolve_col(rows_df, logical, "_b")
            out[logical] = rows_df[actual].values if actual else pd.NA
        return out

    strict_joined_df = _build_joined_df(merged.loc[both_mask])
    relaxed_joined_df = _build_joined_df(relaxed_both) if not relaxed_both.empty else pd.DataFrame(columns=strict_joined_df.columns)
    relaxed_del_joined_df = _build_joined_df(relaxed_del_both) if not relaxed_del_both.empty else pd.DataFrame(columns=strict_joined_df.columns)
    relaxed_startoff_joined_df = _build_joined_df(relaxed_startoff_both) if not relaxed_startoff_both.empty else pd.DataFrame(columns=strict_joined_df.columns)
    joined_out_df = pd.concat([strict_joined_df, relaxed_joined_df, relaxed_del_joined_df, relaxed_startoff_joined_df], ignore_index=True)

    # Apply the canonical column order from the spec
    ordered = [c for c in JOINED_OUTPUT_COLS if c in joined_out_df.columns]
    extra = [c for c in joined_out_df.columns if c not in ordered]
    joined_out_df = joined_out_df[ordered + extra]

    out_joined = os.path.join(args.output_dir, "compare_joined.tsv")
    joined_out_df.to_csv(out_joined, sep="\t", index=False)
    print(f"Wrote joined rows ({len(joined_out_df):,} rows) → {out_joined}")


if __name__ == "__main__":
    main()
