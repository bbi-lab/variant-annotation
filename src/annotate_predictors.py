"""Annotate variants with pre-computed in-silico pathogenicity scores from computational predictors.

Supported scores (hg38 / GRCh38 only):

  REVEL (Rare Exome Variant Ensemble Learner)
  --revel-file revel_hg38.tsv.gz
    Source: https://sites.google.com/site/revelgenomics/downloads
    Range:  0–1 (higher = more likely pathogenic)
    Scope:  missense SNVs only
    --revel-mode {coordinate,transcript,aa} (default: transcript) controls how
    REVEL's one-row-per-overlapping-transcript layout is disambiguated at a
    given genomic position:
      coordinate: take the maximum score across every matching row,
        regardless of transcript. Since a genomic position can be missense
        on one transcript and non-missense (e.g. synonymous) on the
        transcript this project annotates, this mode can attach a REVEL
        score to a variant this project itself calls non-missense.
      transcript (default): additionally require the row's Ensembl
        transcript ID to match this variant's own transcript (from
        --mapped-hgvs-c-transcript-col, mapped RefSeq → Ensembl via
        --revel-mane-file or --revel-transcript-mapping-file). Requires the
        extended REVEL file layout (see "Data file preparation" below).
      aa: additionally require the row's aaref/aaalt columns to equal this
        variant's own amino-acid ref/alt (from --mapped-hgvs-p-ref-col /
        --mapped-hgvs-p-alt-col). Since REVEL's rows are always missense
        (aaref != aaalt), a variant with aa_ref == aa_alt (synonymous) can
        never match — this reproduces the join used by the original
        pipeline notebook. Also requires the extended REVEL file layout.

  AlphaMissense (Google DeepMind)
  --alphamissense-file AlphaMissense_hg38.tsv.gz
    Source: https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
    Index:  generated locally with ``tabix -s 1 -b 2 -e 2 AlphaMissense_hg38.tsv.gz``
    Range:  0–1 (higher = more likely pathogenic)
    Scope:  missense SNVs only

  MutPred2 (via MP2 properties file — preferred)
  --mutpred2-properties-file data_frame_missense_variants_MP2_properties.csv.gz
    A MaveDB-derived per-DNA-variant properties table (columns include
    mavedb_variant_urn, Chrom, hg38_start, hg38_end, ref_allele, alt_allele,
    gene_symbol, AA, "MutPred2 score"). --mutpred2-properties-join-key
    selects how it's joined:
      genomic (default): variant_urn + genomic coordinates, per reverse-
        translation candidate, so scores are pipe-aligned (see below).
      gene-aa: gene symbol + amino-acid substitution (e.g. "T2A", built from
        --gene-symbol-col and the --mapped-hgvs-p-ref/-start/-alt-col
        columns), matching MutPred2's protein-level scope. Evaluated once per
        row, then the same score is duplicated across every pipe slot, since
        the join doesn't vary per DNA candidate. --mutpred2-gene-aa-long-indels
        controls whether candidates with a genomic ref/alt allele over 3bp
        (from --mapped-hgvs-g-ref-col/-alt-col) still get that score, or an
        empty slot instead. --mutpred2-gene-symbol-map-file optionally remaps
        input gene symbols that don't match the properties file's naming
        (e.g. "CALM1_2_3" -> "CALM1") before the lookup.

  MutPred2 (via dbNSFP — legacy fallback)
  --dbnsfp-file dbNSFP5.3.1a_grch38.gz
    Source: https://sites.google.com/site/jpopgen/dbNSFP
    The file and its .tbi index are available pre-built (see run script for URLs).
    Range:  0–1 (higher = more likely pathogenic)
    Scope:  missense SNVs only; returns max score across transcripts.
    Ignored if --mutpred2-properties-file is also given.

At least one of --revel-file, --alphamissense-file, --mutpred2-properties-file,
--dbnsfp-file, --revel-training-file, or --mutpred2-training-file is required.


Performance
-----------

Tabix lookups (REVEL, AlphaMissense, dbNSFP) use pysam's ``TabixFile`` when
pysam is installed, avoiding a subprocess spawn per lookup; a subprocess
fallback runs otherwise. --max-workers N (N > 1) runs those lookups
concurrently across threads within this one process: a first pass over the
input collects all unique SNVs, looks them up in a thread pool, then a
second pass streams the input and annotates from the now-warm cache. This
gets the throughput benefit of running multiple processes over split input
without the split/merge step or the repeated cache-loading cost.


Training-set overlap (optional)
--------------------------------

  --revel-training-file data/revel_training_variants.tsv
    List of REVEL training variants. Joined on hg38 genomic coordinates
    (chromosome, hg38_start, hg38_end, ref_allele, alt_allele) against the
    --mapped-hgvs-g-*-col columns. Produces revel.train, pipe-aligned per
    DNA reverse-translation candidate like revel.score.

  --mutpred2-training-file data/mutpred2_training_variants.tsv
    List of MutPred2 training variants. MutPred2 is a protein-level model,
    so the join key is gene symbol + unqualified protein HGVS rather than
    genomic coordinates: either an 'hgvs_p' column (transcript/protein-
    accession-qualified, e.g. 'NP_000546.1:p.Asn1643His' — gene symbols are
    then ignored), or 'gene_symbol' + 'unqualified_hgvs_p' columns, matched
    against --gene-symbol-col and the part of --mapped-hgvs-p-col following
    the colon. Produces mutpred2.train, duplicated across all DNA candidates
    in the row since the match doesn't vary per candidate.


Data file preparation
---------------------

AlphaMissense is already bgzipped; generate the tabix index locally::

    tabix -s 1 -b 2 -e 2 AlphaMissense_hg38.tsv.gz

dbNSFP GRCh38 variant file and its .tbi index are available pre-built
(download both the .gz and .tbi).

For REVEL, download revel_with_transcript_ids.csv.zip from the link above,
unzip it, then run::

    tail -n +2 revel_with_transcript_ids.csv \\
      | awk -F',' 'NF>=9 && $3!="" && $3!="." \\
                   {print $1"\\t"$3"\\t"$4"\\t"$5"\\t"$8"\\t"$6"\\t"$7"\\t"$9}' \\
      | (printf '#chr\\tpos\\tref\\talt\\trevel_score\\taaref\\taaalt\\tensembl_transcriptid\\n'; sort -k1,1V -k2,2n) \\
      | bgzip > revel_hg38.tsv.gz
    tabix -s 1 -b 2 -e 2 -S 1 revel_hg38.tsv.gz

This extended 8-column layout keeps chr/pos/ref/alt/revel_score in the same
first five columns as before (so --revel-mode coordinate needs nothing
else), and appends aaref/aaalt/ensembl_transcriptid, which --revel-mode
transcript and --revel-mode aa require. A revel_hg38.tsv.gz produced by the
older 5-column command (chr/pos/ref/alt/revel_score only) still works with
--revel-mode coordinate.

--revel-mode transcript additionally needs a RefSeq → Ensembl transcript
mapping, since REVEL's file identifies transcripts by Ensembl ID while this
project's own transcript column is RefSeq. Pass either:

  --revel-mane-file MANE.GRCh38.*.summary.txt.gz
    NCBI MANE summary file (https://ftp.ncbi.nlm.nih.gov/refseq/MANE/MANE_human/current/).
    Restricted to MANE Select / MANE Plus Clinical transcripts, which are
    sequence-identical between RefSeq and Ensembl. Same file format as
    src/remap_transcript_ids.py's --mane-file.
  --revel-transcript-mapping-file FILE
    Custom two-column TSV/CSV with headers source_id, target_id (RefSeq → Ensembl).

The resulting file has five tab-separated columns::

    #chr  pos  ref  alt  revel_score


Output columns
--------------

  revel.score                — REVEL score string (empty for non-SNV / non-missense)
  alphamissense.pathogenicity — AlphaMissense pathogenicity score (0–1)
  alphamissense.class        — likely_benign / ambiguous / likely_pathogenic
  mutpred2.score             — MutPred2 score
  revel.train                — "true"/"false": in the REVEL training set (--revel-training-file)
  mutpred2.train             — "true"/"false": in the MutPred2 training set (--mutpred2-training-file)

For rows with pipe-delimited genomic HGVS values the REVEL and AlphaMissense
columns are pipe-aligned to match the input candidate positions.

mutpred2.score's shape depends on the source:
  --mutpred2-properties-file (--mutpred2-properties-join-key genomic,
    the default): looked up per reverse-translation candidate (variant_urn +
    genomic coordinates), so it is pipe-aligned like REVEL.
  --mutpred2-properties-file (--mutpred2-properties-join-key gene-aa): a
    single score computed once for the row (gene + AA-substitution doesn't
    vary across DNA reverse-translation candidates), then duplicated across
    every pipe slot, so it is pipe-aligned like the genomic join.
  --dbnsfp-file (legacy): MutPred2 is treated as a protein-level model there,
    so all reverse-translation candidates encode the same amino acid
    substitution and a single score (the maximum across candidates) is
    emitted, with no pipes.
"""

from __future__ import annotations

import argparse
import csv
import logging
import os
import re
import subprocess
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from itertools import islice
from pathlib import Path
from typing import Any, Optional

import pandas as pd

from src.remap_transcript_ids import (
    _build_base_mapping,
    _build_custom_mapping,
    _build_mane_mapping,
)

logger = logging.getLogger(__name__)

# Thread-local pysam.TabixFile handles; avoids spawning a subprocess per lookup.
_tabix_local = threading.local()

NC_TO_CHROM_GRCH38: dict[str, str] = {
    "NC_000001.11": "1",
    "NC_000002.12": "2",
    "NC_000003.12": "3",
    "NC_000004.12": "4",
    "NC_000005.10": "5",
    "NC_000006.12": "6",
    "NC_000007.14": "7",
    "NC_000008.11": "8",
    "NC_000009.12": "9",
    "NC_000010.11": "10",
    "NC_000011.10": "11",
    "NC_000012.12": "12",
    "NC_000013.11": "13",
    "NC_000014.9": "14",
    "NC_000015.10": "15",
    "NC_000016.10": "16",
    "NC_000017.11": "17",
    "NC_000018.12": "18",
    "NC_000019.10": "19",
    "NC_000020.11": "20",
    "NC_000021.9": "21",
    "NC_000022.11": "22",
    "NC_000023.11": "X",
    "NC_000024.10": "Y",
    "NC_012920.1": "MT",
}

REVEL_COLS = ["revel.score"]
ALPHAMISSENSE_COLS = ["alphamissense.pathogenicity", "alphamissense.class"]
DBNSFP_COLS = ["mutpred2.score"]
REVEL_TRAIN_COLS = ["revel.train"]
MUTPRED2_TRAIN_COLS = ["mutpred2.train"]

# Module-level cache so the dbNSFP header is read at most once per file path.
_dbnsfp_col_index_cache: dict[str, dict[str, int]] = {}


# ---------------------------------------------------------------------------
# Low-level helpers
# ---------------------------------------------------------------------------

def _chrom_candidates(chrom: str) -> list[str]:
    """Return both "1" and "chr1" variants so we handle either file convention."""
    candidates = [chrom]
    if chrom.startswith("chr"):
        candidates.append(chrom[3:])
    else:
        candidates.append(f"chr{chrom}")
    return list(dict.fromkeys(candidates))


def _get_tabix_handle(path: Path) -> Optional[Any]:
    """Return a thread-local ``pysam.TabixFile`` handle, opening it lazily.

    Returns ``None`` when pysam is unavailable; the caller falls back to
    spawning a subprocess in that case.
    """
    try:
        import pysam  # type: ignore[import]
    except ImportError:
        return None
    handles: Optional[dict[str, Any]] = getattr(_tabix_local, "handles", None)
    if handles is None:
        _tabix_local.handles = {}  # type: ignore[attr-defined]
        handles = _tabix_local.handles
    key = str(path.resolve())
    if key not in handles:
        handles[key] = pysam.TabixFile(str(path))
    return handles[key]


def _run_tabix(path: Path, chrom: str, pos: int) -> list[str]:
    """Return non-comment lines from a tabix point query at *chrom*:*pos*.

    Uses a thread-local ``pysam.TabixFile`` handle when pysam is available
    (avoids spawning a ``tabix`` subprocess per lookup, and lets concurrent
    lookups run in parallel across threads within one process); falls back
    to shelling out to the ``tabix`` binary otherwise.
    """
    handle = _get_tabix_handle(path)
    if handle is not None:
        try:
            return list(handle.fetch(chrom, pos - 1, pos))
        except (ValueError, KeyError):
            return []
    region = f"{chrom}:{pos}-{pos}"
    proc = subprocess.run(
        ["tabix", str(path), region],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        return []
    return [line for line in proc.stdout.splitlines() if line and not line.startswith("#")]


def _snv_from_hgvs_g(
    hgvs_g: str,
    nc_to_chrom: dict[str, str],
) -> Optional[tuple[str, int, str, str]]:
    """Parse a single-base SNV from a genomic HGVS string.

    Returns ``(chrom, pos, ref, alt)`` for recognised SNVs, or ``None`` for
    indels, multi-base substitutions, or unknown NC accessions.
    """
    m = re.match(
        r"^(NC_\d+\.\d+):g\.(\d+)([ACGTacgt]+)>([ACGTacgt]+)$",
        hgvs_g.strip(),
    )
    if not m:
        return None
    chrom = nc_to_chrom.get(m.group(1))
    if not chrom:
        return None
    ref = m.group(3).upper()
    alt = m.group(4).upper()
    if len(ref) != 1 or len(alt) != 1:
        return None  # multi-base MNV, not a true SNV
    return chrom, int(m.group(2)), ref, alt


def _get_dbnsfp_col_indices(path: Path) -> dict[str, int]:
    """Return column-name → 0-based-index mapping from the dbNSFP header.

    Uses ``tabix -H`` to retrieve the header without scanning the entire file.
    Result is cached on *path* so subsequent calls are instantaneous.
    """
    key = str(path)
    if key in _dbnsfp_col_index_cache:
        return _dbnsfp_col_index_cache[key]

    proc = subprocess.run(
        ["tabix", "-H", str(path)],
        capture_output=True,
        text=True,
        check=False,
    )
    header_line: Optional[str] = None
    for line in proc.stdout.splitlines():
        if line.startswith("#"):
            header_line = line
            break

    if header_line is None:
        raise ValueError(
            f"No header line found in {path} via 'tabix -H'. "
            "Ensure the file is tabix-indexed and has a '#'-prefixed header."
        )

    col_names = header_line.lstrip("#").split("\t")
    result = {name.strip(): i for i, name in enumerate(col_names)}
    _dbnsfp_col_index_cache[key] = result
    return result


def _split_pipe(value: str) -> list[str]:
    raw = value or ""
    if "|" not in raw:
        return [raw.strip()]
    return [part.strip() for part in raw.split("|")]


AA_SUBSTITUTION_RE = re.compile(r"^[A-Za-z]$")


def _build_aa_substitution(ref: str, pos: str, alt: str) -> Optional[str]:
    """Build a short-form AA substitution string (e.g. ``T2A``) from split columns.

    Returns ``None`` unless *ref* and *alt* are each a single amino-acid
    letter, *pos* is non-empty, and *ref* != *alt* — matching the properties
    file's missense-only ``AA`` column and naturally excluding stop/
    frameshift/indel/synonymous rows. Mirrors the convention used in
    ``src/find_mp2_unannotated.py``.
    """
    ref = ref.strip().upper()
    alt = alt.strip().upper()
    pos = pos.strip()
    if not (AA_SUBSTITUTION_RE.match(ref) and AA_SUBSTITUTION_RE.match(alt) and pos) or ref == alt:
        return None
    return f"{ref}{pos}{alt}"


def _unqualify_hgvs_p(value: str) -> str:
    """Strip a leading transcript/protein accession from a p. HGVS string.

    ``NP_000546.1:p.Asn1643His`` -> ``p.Asn1643His``. Values with no colon
    (already unqualified) are returned stripped and unchanged.
    """
    value = value.strip()
    if ":" in value:
        return value.rsplit(":", 1)[1].strip()
    return value


# ---------------------------------------------------------------------------
# Per-tool lookups
# ---------------------------------------------------------------------------


def _lookup_mutpred2(
    path: Path,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    cache: dict[tuple[str, int, str, str], Optional[str]],
) -> Optional[str]:
    """Return the maximum MutPred2 score string for an SNV from dbNSFP, or ``None``.

    dbNSFP stores multiple scores per allele (one per protein/transcript),
    semicolon-separated.  A period ``"."`` indicates no score for that entry.
    We return the maximum non-null score across all entries.

    Expected dbNSFP column names (resolved dynamically from the file header):
      chr, pos(1-based), ref, alt, MutPred2_score
    """
    key = (chrom, pos, ref, alt)
    if key in cache:
        return cache[key]

    try:
        col_idx = _get_dbnsfp_col_indices(path)
    except ValueError as exc:
        logger.warning("Cannot read dbNSFP header: %s", exc)
        cache[key] = None
        return None

    chr_col = col_idx.get("chr")
    pos_col = col_idx.get("pos(1-based)")
    ref_col = col_idx.get("ref")
    alt_col = col_idx.get("alt")
    score_col = col_idx.get("MutPred2_score")

    missing = [
        name
        for name, idx in [
            ("chr", chr_col),
            ("pos(1-based)", pos_col),
            ("ref", ref_col),
            ("alt", alt_col),
            ("MutPred2_score", score_col),
        ]
        if idx is None
    ]
    if missing:
        logger.warning(
            "dbNSFP column(s) not found: %s. "
            "Verify that --dbnsfp-file is a dbNSFP GRCh38 variant file.",
            ", ".join(missing),
        )
        cache[key] = None
        return None

    best: Optional[float] = None
    for chrom_try in _chrom_candidates(chrom):
        lines = _run_tabix(path, chrom_try, pos)
        for line in lines:
            fields = line.split("\t")
            try:
                if int(fields[pos_col]) != pos:  # type: ignore[index]
                    continue
            except (ValueError, IndexError):
                continue
            if fields[ref_col].upper() != ref or fields[alt_col].upper() != alt:  # type: ignore[index]
                continue
            if score_col >= len(fields):  # type: ignore[operator]
                continue
            for part in fields[score_col].split(";"):  # type: ignore[index]
                part = part.strip()
                if part in (".", "", "NA"):
                    continue
                try:
                    score = float(part)
                except ValueError:
                    continue
                if best is None or score > best:
                    best = score
        if best is not None:
            break

    result = f"{best:.4f}" if best is not None else None
    cache[key] = result
    return result


def _lookup_mutpred2_from_properties_file(
    file_cache: dict[tuple[str, str, int, int, str, str], str],
    variant_urn: str,
    chrom: str,
    start: str,
    stop: str,
    ref: str,
    alt: str,
) -> Optional[str]:
    """Return the MutPred2 score for one DNA reverse-translation candidate.

    Keyed on ``(variant_urn, chrom, start, stop, ref, alt)`` against a cache
    built by :func:`_load_mutpred2_properties_file_cache`. *chrom* is tried
    both with and without a ``chr`` prefix since the input file and the CVFG
    genomic columns may not agree on convention.
    """
    if not (variant_urn and chrom and start and stop and ref and alt):
        return None
    try:
        start_i = int(start)
        stop_i = int(stop)
    except ValueError:
        return None

    ref_u = ref.strip().upper()
    alt_u = alt.strip().upper()
    for chrom_try in _chrom_candidates(chrom.strip()):
        score = file_cache.get((variant_urn, chrom_try, start_i, stop_i, ref_u, alt_u))
        if score is not None:
            return score
    return None


def _lookup_mutpred2_from_gene_aa_cache(
    cache: dict[tuple[str, str], str],
    gene: str,
    aa: Optional[str],
) -> Optional[str]:
    """Return the MutPred2 score for one (gene, AA-substitution) key, or ``None``.

    Alternative to :func:`_lookup_mutpred2_from_properties_file` for
    ``--mutpred2-properties-join-key gene-aa``. Evaluated once per row rather
    than per DNA reverse-translation candidate, since the join key doesn't
    vary across candidates.
    """
    if not (gene and aa):
        return None
    return cache.get((gene.strip(), aa))


def _lookup_revel_train(
    training_set: set[tuple[str, int, int, str, str]],
    chrom: str,
    start: str,
    stop: str,
    ref: str,
    alt: str,
) -> bool:
    """Return whether one DNA reverse-translation candidate is a REVEL training variant.

    Keyed on genomic coordinates ``(chrom, start, stop, ref, alt)`` against a
    set built by :func:`_load_revel_training_variants`. *chrom* is tried both
    with and without a ``chr`` prefix, as elsewhere in this module.
    """
    if not (chrom and start and stop and ref and alt):
        return False
    try:
        start_i = int(float(start))
        stop_i = int(float(stop))
    except ValueError:
        return False

    ref_u = ref.strip().upper()
    alt_u = alt.strip().upper()
    return any(
        (chrom_try, start_i, stop_i, ref_u, alt_u) in training_set
        for chrom_try in _chrom_candidates(chrom.strip())
    )


def _lookup_mutpred2_train(
    schema: str,
    training_set: set,
    gene_symbol: str,
    mapped_hgvs_p: str,
) -> bool:
    """Return whether a protein variant is a MutPred2 training variant.

    MutPred2 is protein-level, so this is evaluated once per row (not per DNA
    candidate). *mapped_hgvs_p* may itself be pipe-delimited (e.g. multiple
    transcript mappings for the same protein change); a match on any segment
    counts as a match for the whole row. When *schema* is ``"qualified"``,
    segments are matched verbatim against ``training_set`` (a set of
    transcript/protein-accession-qualified hgvs_p strings) and *gene_symbol*
    is ignored. When *schema* is ``"unqualified"``, segments are stripped of
    any accession prefix and matched together with *gene_symbol* against
    ``training_set`` (a set of ``(gene_symbol, unqualified_hgvs_p)`` tuples).
    """
    candidates = [c for c in _split_pipe(mapped_hgvs_p) if c]
    if not candidates:
        return False
    if schema == "qualified":
        return any(c in training_set for c in candidates)
    gene = gene_symbol.strip()
    if not gene:
        return False
    return any((gene, _unqualify_hgvs_p(c)) in training_set for c in candidates)


def _transcript_ids_match(revel_transcript_id: str, target_transcript_id: str) -> bool:
    """Compare an Ensembl transcript accession from REVEL against the resolved target.

    Tries an exact (versioned) match first, then falls back to the
    version-stripped base accession, since REVEL and a MANE mapping may pin
    different versions of the same underlying transcript.
    """
    if not revel_transcript_id or not target_transcript_id:
        return False
    if revel_transcript_id == target_transcript_id:
        return True
    base_a = revel_transcript_id.rsplit(".", 1)[0] if "." in revel_transcript_id else revel_transcript_id
    base_b = target_transcript_id.rsplit(".", 1)[0] if "." in target_transcript_id else target_transcript_id
    return base_a == base_b


def _resolve_ensembl_transcript(
    refseq_transcript: str,
    mapping: dict[str, str],
    base_mapping: dict[str, str],
) -> Optional[str]:
    """Resolve a RefSeq transcript accession to its Ensembl equivalent, or ``None``.

    Tries the full versioned accession first, then the version-stripped base
    accession. Unlike :func:`src.remap_transcript_ids._remap_accession`, an
    accession absent from both mappings resolves to ``None`` rather than
    being returned unchanged — a REVEL transcript-mode lookup with no known
    Ensembl equivalent must not silently fall back to matching nothing as if
    it were a literal accession.
    """
    if not refseq_transcript:
        return None
    if refseq_transcript in mapping:
        return mapping[refseq_transcript]
    base = refseq_transcript.rsplit(".", 1)[0] if "." in refseq_transcript else refseq_transcript
    return base_mapping.get(base)


_REVEL_MODE_MIN_COLS = {"coordinate": 5, "aa": 7, "transcript": 8}


def _get_revel_header_col_count(path: Path) -> Optional[int]:
    """Return the number of tab-separated columns in *path*'s '#'-prefixed header, or ``None``."""
    proc = subprocess.run(["tabix", "-H", str(path)], capture_output=True, text=True, check=False)
    header_line = next((line for line in proc.stdout.splitlines() if line.startswith("#")), None)
    if header_line is None:
        return None
    return len(header_line.lstrip("#").split("\t"))


def _validate_revel_file_for_mode(path: Path, mode: str) -> None:
    """Raise ``ValueError`` if *path* lacks the columns required by *mode*."""
    required = _REVEL_MODE_MIN_COLS[mode]
    if required <= 5:
        return
    n_cols = _get_revel_header_col_count(path)
    if n_cols is None or n_cols < required:
        extra = "aaref/aaalt/ensembl_transcriptid" if mode == "transcript" else "aaref/aaalt"
        raise ValueError(
            f"REVEL file {path} does not have the extended column layout required by "
            f"--revel-mode {mode} (found {n_cols if n_cols is not None else 'no'} header "
            f"column(s) via 'tabix -H', need at least {required}, including {extra}). "
            "Regenerate it per the REVEL section of docs/annotate_predictors.md, or use "
            "--revel-mode coordinate with the existing file."
        )


def _lookup_revel(
    path: Path,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    cache: dict[tuple, Optional[str]],
    *,
    mode: str = "coordinate",
    transcript_id: Optional[str] = None,
    aa_ref: Optional[str] = None,
    aa_alt: Optional[str] = None,
) -> Optional[str]:
    """Return the maximum REVEL score string for an SNV, or ``None`` if absent.

    REVEL may have multiple rows per position (one per overlapping
    transcript). *mode* controls how those rows are disambiguated:

      "coordinate" (default here; the annotate_predictors CLI default is
        "transcript" — see --revel-mode): take the maximum score across
        every matching row regardless of transcript. Since a genomic
        position can be missense on one transcript and non-missense (e.g.
        synonymous) on another, this mode can attach a score to a variant
        whose own annotated consequence is non-missense.
      "transcript": additionally require the row's Ensembl transcript-ID
        column (field 7) to match *transcript_id* (see
        :func:`_transcript_ids_match`).
      "aa": additionally require the row's aaref/aaalt columns (fields 5, 6)
        to equal *aa_ref*/*aa_alt* (case-insensitive). Since REVEL's rows
        are always missense (aaref != aaalt), a variant with aa_ref ==
        aa_alt (synonymous) can never match.

    Expected prepared-file column layout (tab-separated, 0-indexed):
      0: chr  1: pos  2: ref  3: alt  4: revel_score
      Extended (required for "transcript"/"aa" modes — see
      docs/annotate_predictors.md for the preparation command):
      5: aaref  6: aaalt  7: ensembl_transcriptid
    """
    aa_ref = aa_ref.upper() if aa_ref else None
    aa_alt = aa_alt.upper() if aa_alt else None
    key = (mode, chrom, pos, ref, alt, transcript_id or "", aa_ref or "", aa_alt or "")
    if key in cache:
        return cache[key]

    if mode == "transcript" and not transcript_id:
        cache[key] = None
        return None
    if mode == "aa" and not (aa_ref and aa_alt):
        cache[key] = None
        return None

    best: Optional[float] = None
    for chrom_try in _chrom_candidates(chrom):
        lines = _run_tabix(path, chrom_try, pos)
        for line in lines:
            fields = line.split("\t")
            if len(fields) < 5:
                continue
            try:
                r_pos = int(fields[1])
            except ValueError:
                continue
            if r_pos != pos:
                continue
            if fields[2].upper() != ref or fields[3].upper() != alt:
                continue
            if mode == "transcript" and (
                len(fields) < 8 or not _transcript_ids_match(fields[7].strip(), transcript_id)  # type: ignore[arg-type]
            ):
                continue
            if mode == "aa" and (
                len(fields) < 7 or fields[5].strip().upper() != aa_ref or fields[6].strip().upper() != aa_alt
            ):
                continue
            try:
                score = float(fields[4])
            except ValueError:
                continue
            if best is None or score > best:
                best = score
        if best is not None:
            break

    result = f"{best:.4f}" if best is not None else None
    cache[key] = result
    return result


def _lookup_alphamissense(
    path: Path,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    cache: dict[tuple[str, int, str, str], Optional[tuple[str, str]]],
) -> Optional[tuple[str, str]]:
    """Return ``(am_pathogenicity, am_class)`` for an SNV, or ``None`` if absent.

    When multiple transcript entries exist for the same allele (which is typical
    in AlphaMissense), the entry with the highest pathogenicity score is returned.

    AlphaMissense column layout (tab-separated, 0-indexed):
      0: CHROM  1: POS  2: REF  3: ALT  4: genome  5: uniprot_id
      6: transcript_id  7: protein_variant  8: am_pathogenicity  9: am_class
    """
    key = (chrom, pos, ref, alt)
    if key in cache:
        return cache[key]

    best_score: Optional[float] = None
    best_class = ""
    for chrom_try in _chrom_candidates(chrom):
        lines = _run_tabix(path, chrom_try, pos)
        for line in lines:
            fields = line.split("\t")
            if len(fields) < 10:
                continue
            try:
                r_pos = int(fields[1])
            except ValueError:
                continue
            if r_pos != pos:
                continue
            if fields[2].upper() != ref or fields[3].upper() != alt:
                continue
            try:
                score = float(fields[8])
            except ValueError:
                continue
            if best_score is None or score > best_score:
                best_score = score
                best_class = fields[9].strip()
        if best_score is not None:
            break

    result: Optional[tuple[str, str]] = (
        (f"{best_score:.4f}", best_class) if best_score is not None else None
    )
    cache[key] = result
    return result


# ---------------------------------------------------------------------------
# File-based caches
# ---------------------------------------------------------------------------


def _load_revel_file_cache(path: str) -> dict[str, str]:
    """Load a two-column TSV (hgvs, revel.score) into an HGVS → score dict."""
    cache: dict[str, str] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            h = row.get("hgvs", "").strip()
            s = row.get("revel.score", "").strip()
            if h and s:
                cache[h] = s
    logger.info("Loaded %d REVEL file-cache entries from %s", len(cache), path)
    return cache


def _load_alphamissense_file_cache(path: str) -> dict[str, tuple[str, str]]:
    """Load a three-column TSV into an HGVS → (pathogenicity, class) dict."""
    cache: dict[str, tuple[str, str]] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            h  = row.get("hgvs", "").strip()
            ap = row.get("alphamissense.pathogenicity", "").strip()
            ac = row.get("alphamissense.class", "").strip()
            if h and (ap or ac):
                cache[h] = (ap, ac)
    logger.info("Loaded %d AlphaMissense file-cache entries from %s", len(cache), path)
    return cache


MUTPRED2_PROPERTIES_REQUIRED_COLS = [
    "mavedb_variant_urn",
    "Chrom",
    "hg38_start",
    "hg38_end",
    "ref_allele",
    "alt_allele",
    "MutPred2 score",
]


def _load_mutpred2_properties_file_cache(
    path: str,
) -> dict[tuple[str, str, int, int, str, str], str]:
    """Load MutPred2 scores from a MaveDB MP2-properties CSV into a lookup dict.

    Expects (at least) the columns ``mavedb_variant_urn``, ``Chrom``,
    ``hg38_start``, ``hg38_end``, ``ref_allele``, ``alt_allele``, and
    ``MutPred2 score``. The file describes one row per DNA variant, so the
    same ``mavedb_variant_urn`` (a MaveDB protein- or DNA-level variant) can
    repeat across multiple reverse-translation candidates; each is keyed
    separately on ``(variant_urn, chrom, start, stop, ref, alt)``.

    Transparently handles gzip-compressed input (``.gz`` suffix). Uses
    pandas' C parser (reading only the required columns) rather than a
    per-row Python loop, since this file commonly has well over a million
    rows and the naive loop dominates startup time.
    """
    try:
        header_cols = pd.read_csv(path, nrows=0).columns
    except pd.errors.EmptyDataError:
        raise ValueError(f"MutPred2 properties file is empty: {path}")

    missing = [c for c in MUTPRED2_PROPERTIES_REQUIRED_COLS if c not in header_cols]
    if missing:
        raise ValueError(
            f"MutPred2 properties file {path} is missing column(s): {', '.join(missing)}"
        )

    df = pd.read_csv(
        path,
        usecols=MUTPRED2_PROPERTIES_REQUIRED_COLS,
        dtype=str,
        na_filter=False,
    )
    n_rows = len(df)

    urn = df["mavedb_variant_urn"].str.strip()
    score = df["MutPred2 score"].str.strip()
    start = pd.to_numeric(df["hg38_start"], errors="coerce")
    stop = pd.to_numeric(df["hg38_end"], errors="coerce")

    valid = (urn != "") & (score != "") & start.notna() & stop.notna()
    n_skipped = n_rows - int(valid.sum())

    chrom = df["Chrom"].str.strip()[valid]
    ref = df["ref_allele"].str.strip().str.upper()[valid]
    alt = df["alt_allele"].str.strip().str.upper()[valid]

    keys = zip(
        urn[valid].tolist(),
        chrom.tolist(),
        start[valid].astype(int).tolist(),
        stop[valid].astype(int).tolist(),
        ref.tolist(),
        alt.tolist(),
    )
    cache: dict[tuple[str, str, int, int, str, str], str] = dict(zip(keys, score[valid].tolist()))

    logger.info(
        "Loaded %d MutPred2 file-cache entries from %s (%d rows scanned, %d skipped)",
        len(cache), path, n_rows, n_skipped,
    )
    return cache


MUTPRED2_PROPERTIES_GENE_AA_REQUIRED_COLS = ["gene_symbol", "AA", "MutPred2 score"]


def _load_mutpred2_properties_gene_aa_cache(path: str) -> dict[tuple[str, str], str]:
    """Load MutPred2 scores from a MaveDB MP2-properties CSV into a (gene, AA) -> score dict.

    Alternative to :func:`_load_mutpred2_properties_file_cache` for
    ``--mutpred2-properties-join-key gene-aa``. MutPred2 is a protein-level
    model, so a variant only needs one score per (gene, amino-acid
    substitution) regardless of which mavedb_variant_urn or DNA reverse-
    translation candidate it came from in the properties file. Keyed on the
    ``gene_symbol`` column (the file's ``Gene`` column is ignored).

    Transparently handles gzip-compressed input (``.gz`` suffix).
    """
    try:
        header_cols = pd.read_csv(path, nrows=0).columns
    except pd.errors.EmptyDataError:
        raise ValueError(f"MutPred2 properties file is empty: {path}")

    missing = [c for c in MUTPRED2_PROPERTIES_GENE_AA_REQUIRED_COLS if c not in header_cols]
    if missing:
        raise ValueError(
            f"MutPred2 properties file {path} is missing column(s): {', '.join(missing)}"
        )

    df = pd.read_csv(
        path,
        usecols=MUTPRED2_PROPERTIES_GENE_AA_REQUIRED_COLS,
        dtype=str,
        na_filter=False,
    )
    n_rows = len(df)

    gene_symbol = df["gene_symbol"].str.strip()
    aa = df["AA"].str.strip()
    score = df["MutPred2 score"].str.strip()

    valid = (gene_symbol != "") & (aa != "") & (score != "")
    n_skipped = n_rows - int(valid.sum())

    keys = zip(gene_symbol[valid].tolist(), aa[valid].tolist())
    cache: dict[tuple[str, str], str] = dict(zip(keys, score[valid].tolist()))

    logger.info(
        "Loaded %d MutPred2 (gene, AA) file-cache entries from %s (%d rows scanned, %d skipped)",
        len(cache), path, n_rows, n_skipped,
    )
    return cache


MUTPRED2_GENE_SYMBOL_MAP_REQUIRED_COLS = ["gene_symbol", "mutpred2_gene_symbol"]


def _load_mutpred2_gene_symbol_map(path: str) -> dict[str, str]:
    """Load a TSV mapping input gene symbols to MutPred2-properties-file gene symbols.

    For ``--mutpred2-properties-join-key gene-aa``, when a gene is named
    differently between the two files (e.g. a combined-target dataset's
    ``CALM1_2_3`` in the input vs. ``CALM1`` in the properties file's
    ``gene_symbol`` column), this remaps the input gene symbol before the
    (gene, AA) cache lookup. Only genes that actually differ need an entry;
    any gene symbol absent from the map is looked up unchanged.

    Expects a header with ``gene_symbol`` (as it appears via
    ``--gene-symbol-col`` in the input) and ``mutpred2_gene_symbol`` (as it
    appears in the properties file). A repeated ``gene_symbol`` keeps its
    last row's mapping.
    """
    mapping: dict[str, str] = {}
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"MutPred2 gene-symbol-map file is empty: {path}")
        missing = [c for c in MUTPRED2_GENE_SYMBOL_MAP_REQUIRED_COLS if c not in reader.fieldnames]
        if missing:
            raise ValueError(
                f"MutPred2 gene-symbol-map file {path} is missing column(s): {', '.join(missing)}"
            )

        for row in reader:
            src = (row.get("gene_symbol") or "").strip()
            dst = (row.get("mutpred2_gene_symbol") or "").strip()
            if src and dst:
                mapping[src] = dst

    logger.info("Loaded %d MutPred2 gene-symbol mapping(s) from %s", len(mapping), path)
    return mapping


REVEL_TRAINING_REQUIRED_COLS = [
    "chromosome",
    "hg38_start",
    "hg38_end",
    "ref_allele",
    "alt_allele",
]


def _load_revel_training_variants(path: str) -> set[tuple[str, int, int, str, str]]:
    """Load a REVEL training-variant TSV into a set of genomic-coordinate keys.

    Expects (at least) the columns ``chromosome``, ``hg38_start``,
    ``hg38_end``, ``ref_allele``, and ``alt_allele``. Any other columns (e.g.
    ``gene_symbol``, ``unqualified_hgvs_p``) are informational only and are
    not used for matching. Positions may be written as floats (e.g.
    ``3476169.0``), which is tolerated.
    """
    training_set: set[tuple[str, int, int, str, str]] = set()
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"REVEL training-variants file is empty: {path}")
        missing = [c for c in REVEL_TRAINING_REQUIRED_COLS if c not in reader.fieldnames]
        if missing:
            raise ValueError(
                f"REVEL training-variants file {path} is missing column(s): {', '.join(missing)}"
            )

        for row in reader:
            chrom = (row.get("chromosome") or "").strip()
            ref = (row.get("ref_allele") or "").strip().upper()
            alt = (row.get("alt_allele") or "").strip().upper()
            try:
                start = int(float(row.get("hg38_start") or ""))
                stop = int(float(row.get("hg38_end") or ""))
            except ValueError:
                continue
            if not (chrom and ref and alt):
                continue
            training_set.add((chrom, start, stop, ref, alt))

    logger.info("Loaded %d REVEL training-variant keys from %s", len(training_set), path)
    return training_set


MUTPRED2_TRAINING_QUALIFIED_COL = "hgvs_p"
MUTPRED2_TRAINING_UNQUALIFIED_COLS = ["gene_symbol", "unqualified_hgvs_p"]


def _load_mutpred2_training_variants(path: str) -> tuple[str, set]:
    """Load a MutPred2 training-variant TSV, auto-detecting its join schema.

    If the file has an ``hgvs_p`` column (a transcript/protein-accession-
    qualified HGVS string, e.g. ``NP_000546.1:p.Asn1643His``), matching is
    done on that full qualified string and gene symbols are ignored.
    Otherwise the file must have ``gene_symbol`` and ``unqualified_hgvs_p``
    columns, and matching is done on that pair.

    Returns ``(schema, training_set)`` where *schema* is ``"qualified"`` or
    ``"unqualified"``, and *training_set* holds either qualified hgvs_p
    strings or ``(gene_symbol, unqualified_hgvs_p)`` tuples, respectively.
    """
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"MutPred2 training-variants file is empty: {path}")
        fieldnames = set(reader.fieldnames)

        if MUTPRED2_TRAINING_QUALIFIED_COL in fieldnames:
            schema = "qualified"
            training_set: set = set()
            for row in reader:
                hgvs_p = (row.get(MUTPRED2_TRAINING_QUALIFIED_COL) or "").strip()
                if hgvs_p:
                    training_set.add(hgvs_p)
        elif fieldnames.issuperset(MUTPRED2_TRAINING_UNQUALIFIED_COLS):
            schema = "unqualified"
            training_set = set()
            for row in reader:
                gene = (row.get("gene_symbol") or "").strip()
                hgvs_p = (row.get("unqualified_hgvs_p") or "").strip()
                if gene and hgvs_p:
                    training_set.add((gene, hgvs_p))
        else:
            raise ValueError(
                f"MutPred2 training-variants file {path} must have either an "
                f"'{MUTPRED2_TRAINING_QUALIFIED_COL}' column or all of "
                f"{MUTPRED2_TRAINING_UNQUALIFIED_COLS}"
            )

    logger.info(
        "Loaded %d MutPred2 training-variant keys from %s (schema: %s)",
        len(training_set), path, schema,
    )
    return schema, training_set


# ---------------------------------------------------------------------------
# Row-level annotation
# ---------------------------------------------------------------------------

def annotate_row(
    row: dict[str, str],
    *,
    nc_to_chrom: dict[str, str],
    mapped_hgvs_g_col: str,
    mapped_hgvs_c_col: Optional[str] = None,
    revel_path: Optional[Path],
    alphamissense_path: Optional[Path],
    dbnsfp_path: Optional[Path] = None,
    revel_cache: dict[tuple, Optional[str]],
    revel_mode: str = "coordinate",
    revel_transcript_mapping: Optional[dict[str, str]] = None,
    revel_transcript_base_mapping: Optional[dict[str, str]] = None,
    mapped_hgvs_c_transcript_col: Optional[str] = None,
    am_cache: dict[tuple[str, int, str, str], Optional[tuple[str, str]]],
    mutpred2_cache: Optional[dict[tuple[str, int, str, str], Optional[str]]] = None,
    revel_file_cache: Optional[dict[str, str]] = None,
    am_file_cache: Optional[dict[str, tuple[str, str]]] = None,
    mutpred2_properties_cache: Optional[dict[tuple[str, str, int, int, str, str], str]] = None,
    mutpred2_gene_aa_cache: Optional[dict[tuple[str, str], str]] = None,
    mutpred2_gene_aa_long_indels: str = "annotate",
    mutpred2_gene_symbol_map: Optional[dict[str, str]] = None,
    variant_urn_col: Optional[str] = None,
    mapped_hgvs_g_chromosome_col: Optional[str] = None,
    mapped_hgvs_g_start_col: Optional[str] = None,
    mapped_hgvs_g_stop_col: Optional[str] = None,
    mapped_hgvs_g_ref_col: Optional[str] = None,
    mapped_hgvs_g_alt_col: Optional[str] = None,
    revel_training_set: Optional[set[tuple[str, int, int, str, str]]] = None,
    mutpred2_training_schema: Optional[str] = None,
    mutpred2_training_set: Optional[set] = None,
    gene_symbol_col: Optional[str] = None,
    mapped_hgvs_p_col: Optional[str] = None,
    mapped_hgvs_p_ref_col: Optional[str] = None,
    mapped_hgvs_p_start_col: Optional[str] = None,
    mapped_hgvs_p_alt_col: Optional[str] = None,
) -> dict[str, str]:
    """Return annotation columns for a single row.

    Output values are pipe-aligned to the candidates.  For REVEL and
    AlphaMissense each candidate is resolved by checking the file cache
    (keyed on the c-string from *mapped_hgvs_c_col*) first, then falling back
    to a tabix lookup via the g-string from *mapped_hgvs_g_col*.  Non-SNV
    candidates that are absent from the file cache produce empty strings.

    *revel_mode* (default "coordinate" here; the CLI default is
    "transcript") is passed through to :func:`_lookup_revel` for every
    tabix-backed REVEL candidate — see that function for what each mode
    means. "transcript" mode resolves this row's own RefSeq transcript (from
    *mapped_hgvs_c_transcript_col*) to Ensembl via *revel_transcript_mapping*
    / *revel_transcript_base_mapping* once per row (the transcript doesn't
    vary across DNA reverse-translation candidates). "aa" mode reads the
    row's own amino-acid ref/alt from *mapped_hgvs_p_ref_col* /
    *mapped_hgvs_p_alt_col*, also once per row.

    MutPred2 has three, mutually-exclusive sources. If *mutpred2_properties_cache*
    is given, each reverse-translation candidate is looked up individually
    (variant_urn + genomic coordinates) and the result is pipe-aligned like
    REVEL/AlphaMissense. If *mutpred2_gene_aa_cache* is given instead, the
    join is gene symbol + amino-acid substitution (built from
    *mapped_hgvs_p_ref_col*/*_start_col*/*_alt_col*), evaluated once for the
    row since MutPred2 is protein-level and the key doesn't vary per DNA
    candidate, then the same score is duplicated across every pipe slot so it
    is still pipe-aligned like REVEL/AlphaMissense. If *mutpred2_gene_symbol_map*
    is given, the row's gene symbol is remapped through it (e.g. ``CALM1_2_3``
    -> ``CALM1``) before the (gene, AA) cache lookup; genes absent from the
    map are looked up unchanged. If *mutpred2_gene_aa_long_indels*
    is ``"ignore"``, any candidate whose genomic ref or alt allele (from
    *mapped_hgvs_g_ref_col*/*_alt_col*) exceeds 3bp gets an empty slot instead
    of the row's score, since such a large indel/delins is a poor proxy for
    the assayed protein substitution; the default, ``"annotate"``, applies the
    row's score to every candidate regardless of allele length. Otherwise, if
    *dbnsfp_path* is given, MutPred2 is treated as a protein-level model and a
    single score (the maximum across candidates) is emitted with no pipes.

    Training-set overlap is optional and independent of the score sources
    above. If *revel_training_set* is given, ``revel.train`` is looked up per
    DNA candidate (genomic coordinates) and pipe-aligned like ``revel.score``.
    If *mutpred2_training_set* is given, ``mutpred2.train`` is evaluated once
    for the row (gene symbol + protein HGVS) and the same "true"/"false"
    value is duplicated across all DNA candidates, since the match doesn't
    vary per candidate.
    """
    g_candidates = _split_pipe((row.get(mapped_hgvs_g_col) or "").strip())
    c_candidates = (
        _split_pipe((row.get(mapped_hgvs_c_col) or "").strip())
        if mapped_hgvs_c_col
        else []
    )

    mutpred2_alt_enabled = mutpred2_properties_cache is not None
    mutpred2_gene_aa_enabled = mutpred2_gene_aa_cache is not None
    mutpred2_gene_aa_skip_long_indels = (
        mutpred2_gene_aa_enabled and mutpred2_gene_aa_long_indels == "ignore"
    )
    revel_train_enabled = revel_training_set is not None
    mutpred2_train_enabled = mutpred2_training_set is not None
    need_geno_candidates = mutpred2_alt_enabled or revel_train_enabled or mutpred2_gene_aa_skip_long_indels

    variant_urn = ""
    if mutpred2_alt_enabled:
        variant_urn = (row.get(variant_urn_col) or "").strip() if variant_urn_col else ""

    chrom_candidates: list[str] = []
    start_candidates: list[str] = []
    stop_candidates: list[str] = []
    ref_candidates: list[str] = []
    alt_candidates: list[str] = []
    if need_geno_candidates:
        if mapped_hgvs_g_chromosome_col:
            chrom_candidates = _split_pipe((row.get(mapped_hgvs_g_chromosome_col) or "").strip())
        if mapped_hgvs_g_start_col:
            start_candidates = _split_pipe((row.get(mapped_hgvs_g_start_col) or "").strip())
        if mapped_hgvs_g_stop_col:
            stop_candidates = _split_pipe((row.get(mapped_hgvs_g_stop_col) or "").strip())
        if mapped_hgvs_g_ref_col:
            ref_candidates = _split_pipe((row.get(mapped_hgvs_g_ref_col) or "").strip())
        if mapped_hgvs_g_alt_col:
            alt_candidates = _split_pipe((row.get(mapped_hgvs_g_alt_col) or "").strip())

    n = max(
        len(g_candidates),
        len(c_candidates),
        len(chrom_candidates),
        len(start_candidates),
        len(stop_candidates),
        len(ref_candidates),
        len(alt_candidates),
    )
    g_candidates += [""] * (n - len(g_candidates))
    c_candidates += [""] * (n - len(c_candidates))
    chrom_candidates += [""] * (n - len(chrom_candidates))
    start_candidates += [""] * (n - len(start_candidates))
    stop_candidates += [""] * (n - len(stop_candidates))
    ref_candidates += [""] * (n - len(ref_candidates))
    alt_candidates += [""] * (n - len(alt_candidates))

    revel_enabled = revel_path is not None or revel_file_cache is not None
    am_enabled = alphamissense_path is not None or am_file_cache is not None

    revel_transcript_id: Optional[str] = None
    revel_aa_ref: Optional[str] = None
    revel_aa_alt: Optional[str] = None
    if revel_enabled and revel_path is not None and revel_mode == "transcript":
        raw_transcript = (
            (row.get(mapped_hgvs_c_transcript_col) or "").strip() if mapped_hgvs_c_transcript_col else ""
        )
        if raw_transcript:
            raw_transcript = _split_pipe(raw_transcript)[0]
        revel_transcript_id = _resolve_ensembl_transcript(
            raw_transcript, revel_transcript_mapping or {}, revel_transcript_base_mapping or {}
        )
    elif revel_enabled and revel_path is not None and revel_mode == "aa":
        raw_aa_ref = (row.get(mapped_hgvs_p_ref_col) or "").strip() if mapped_hgvs_p_ref_col else ""
        raw_aa_alt = (row.get(mapped_hgvs_p_alt_col) or "").strip() if mapped_hgvs_p_alt_col else ""
        revel_aa_ref = raw_aa_ref.upper() or None
        revel_aa_alt = raw_aa_alt.upper() or None

    revel_vals: list[str] = []
    am_path_vals: list[str] = []
    am_class_vals: list[str] = []
    mutpred2_vals: list[str] = []
    mutpred2_alt_vals: list[str] = []
    revel_train_vals: list[str] = []

    _mp2_cache: dict[tuple[str, int, str, str], Optional[str]] = (
        mutpred2_cache if mutpred2_cache is not None else {}
    )

    for i in range(n):
        hgvs_g = g_candidates[i]
        hgvs_c = c_candidates[i]
        snv = _snv_from_hgvs_g(hgvs_g, nc_to_chrom) if hgvs_g else None

        if revel_enabled:
            r: Optional[str] = None
            if hgvs_c and revel_file_cache is not None:
                r = revel_file_cache.get(hgvs_c)
            if r is None and snv is not None and revel_path is not None:
                chrom, pos, ref, alt = snv
                r = _lookup_revel(
                    revel_path, chrom, pos, ref, alt, revel_cache,
                    mode=revel_mode,
                    transcript_id=revel_transcript_id,
                    aa_ref=revel_aa_ref,
                    aa_alt=revel_aa_alt,
                )
            revel_vals.append(r or "")

        if am_enabled:
            a: Optional[tuple[str, str]] = None
            if hgvs_c and am_file_cache is not None:
                a = am_file_cache.get(hgvs_c)
            if a is None and snv is not None and alphamissense_path is not None:
                chrom, pos, ref, alt = snv
                a = _lookup_alphamissense(alphamissense_path, chrom, pos, ref, alt, am_cache)
            if a is not None:
                am_path_vals.append(a[0])
                am_class_vals.append(a[1])
            else:
                am_path_vals.append("")
                am_class_vals.append("")

        if mutpred2_alt_enabled:
            m_alt = _lookup_mutpred2_from_properties_file(
                mutpred2_properties_cache,  # type: ignore[arg-type]
                variant_urn,
                chrom_candidates[i],
                start_candidates[i],
                stop_candidates[i],
                ref_candidates[i],
                alt_candidates[i],
            )
            mutpred2_alt_vals.append(m_alt or "")
        elif snv is not None and dbnsfp_path is not None:
            chrom, pos, ref, alt = snv
            m = _lookup_mutpred2(dbnsfp_path, chrom, pos, ref, alt, _mp2_cache)
            if m is not None:
                mutpred2_vals.append(m)

        if revel_train_enabled:
            is_train = _lookup_revel_train(
                revel_training_set,  # type: ignore[arg-type]
                chrom_candidates[i],
                start_candidates[i],
                stop_candidates[i],
                ref_candidates[i],
                alt_candidates[i],
            )
            revel_train_vals.append("true" if is_train else "false")

    sep = "|" if n > 1 else ""
    out: dict[str, str] = {}
    if revel_enabled:
        out["revel.score"] = sep.join(revel_vals)
    if am_enabled:
        out["alphamissense.pathogenicity"] = sep.join(am_path_vals)
        out["alphamissense.class"] = sep.join(am_class_vals)
    if mutpred2_alt_enabled:
        out["mutpred2.score"] = sep.join(mutpred2_alt_vals)
    elif mutpred2_gene_aa_enabled:
        gene = (row.get(gene_symbol_col) or "").strip() if gene_symbol_col else ""
        if mutpred2_gene_symbol_map:
            gene = mutpred2_gene_symbol_map.get(gene, gene)
        ref = (row.get(mapped_hgvs_p_ref_col) or "") if mapped_hgvs_p_ref_col else ""
        pos = (row.get(mapped_hgvs_p_start_col) or "") if mapped_hgvs_p_start_col else ""
        alt = (row.get(mapped_hgvs_p_alt_col) or "") if mapped_hgvs_p_alt_col else ""
        aa = _build_aa_substitution(ref, pos, alt)
        gene_aa_score = _lookup_mutpred2_from_gene_aa_cache(
            mutpred2_gene_aa_cache,  # type: ignore[arg-type]
            gene,
            aa,
        ) or ""
        if mutpred2_gene_aa_skip_long_indels:
            gene_aa_vals = [
                ""
                if max(len(ref_candidates[i]), len(alt_candidates[i])) > 3
                else gene_aa_score
                for i in range(n)
            ]
        else:
            gene_aa_vals = [gene_aa_score] * n
        out["mutpred2.score"] = sep.join(gene_aa_vals)
    elif dbnsfp_path is not None:
        if mutpred2_vals:
            best_mp2 = max(mutpred2_vals, key=float)
        else:
            best_mp2 = ""
        out["mutpred2.score"] = best_mp2
    if revel_train_enabled:
        out["revel.train"] = sep.join(revel_train_vals)
    if mutpred2_train_enabled:
        gene_symbol = (row.get(gene_symbol_col) or "") if gene_symbol_col else ""
        mapped_hgvs_p = (row.get(mapped_hgvs_p_col) or "") if mapped_hgvs_p_col else ""
        is_train = _lookup_mutpred2_train(
            mutpred2_training_schema,  # type: ignore[arg-type]
            mutpred2_training_set,  # type: ignore[arg-type]
            gene_symbol,
            mapped_hgvs_p,
        )
        out["mutpred2.train"] = sep.join(["true" if is_train else "false"] * n)
    return out


# ---------------------------------------------------------------------------
# Parallel tabix prefetch (--max-workers > 1)
# ---------------------------------------------------------------------------

def _collect_needed_snv_keys(
    rows: Any,
    *,
    nc_to_chrom: dict[str, str],
    mapped_hgvs_g_col: str,
    mapped_hgvs_c_col: Optional[str],
    revel_enabled: bool,
    revel_file_cache: Optional[dict[str, str]],
    am_enabled: bool,
    am_file_cache: Optional[dict[str, tuple[str, str]]],
    mutpred2_dbnsfp_enabled: bool,
    revel_mode: str = "coordinate",
    revel_transcript_mapping: Optional[dict[str, str]] = None,
    revel_transcript_base_mapping: Optional[dict[str, str]] = None,
    mapped_hgvs_c_transcript_col: Optional[str] = None,
    mapped_hgvs_p_ref_col: Optional[str] = None,
    mapped_hgvs_p_alt_col: Optional[str] = None,
) -> tuple[set[tuple], set[tuple[str, int, str, str]], set[tuple[str, int, str, str]]]:
    """Scan *rows* and return the unique SNV keys that will need a tabix lookup.

    Mirrors the candidate-splitting logic in :func:`annotate_row`, but only
    for the g./c. HGVS columns that feed REVEL/AlphaMissense/dbNSFP — the
    other lookups (MutPred2-properties-file, training-set overlap) are plain
    dict/set lookups and don't benefit from tabix prefetching.

    *revel_needed* entries are 8-tuples matching :func:`_lookup_revel`'s
    cache key (``mode, chrom, pos, ref, alt, transcript_id, aa_ref, aa_alt``)
    so prefetched results land in the same cache slots the main pass reads.
    """
    revel_needed: set[tuple] = set()
    am_needed: set[tuple[str, int, str, str]] = set()
    mp2_needed: set[tuple[str, int, str, str]] = set()

    for row in rows:
        g_candidates = _split_pipe((row.get(mapped_hgvs_g_col) or "").strip())
        c_candidates = (
            _split_pipe((row.get(mapped_hgvs_c_col) or "").strip())
            if mapped_hgvs_c_col
            else []
        )

        revel_transcript_id: Optional[str] = None
        revel_aa_ref: Optional[str] = None
        revel_aa_alt: Optional[str] = None
        if revel_enabled and revel_mode == "transcript":
            raw_transcript = (
                (row.get(mapped_hgvs_c_transcript_col) or "").strip() if mapped_hgvs_c_transcript_col else ""
            )
            if raw_transcript:
                raw_transcript = _split_pipe(raw_transcript)[0]
            revel_transcript_id = _resolve_ensembl_transcript(
                raw_transcript, revel_transcript_mapping or {}, revel_transcript_base_mapping or {}
            )
        elif revel_enabled and revel_mode == "aa":
            raw_aa_ref = (row.get(mapped_hgvs_p_ref_col) or "").strip() if mapped_hgvs_p_ref_col else ""
            raw_aa_alt = (row.get(mapped_hgvs_p_alt_col) or "").strip() if mapped_hgvs_p_alt_col else ""
            revel_aa_ref = raw_aa_ref.upper() or None
            revel_aa_alt = raw_aa_alt.upper() or None

        for i, hgvs_g in enumerate(g_candidates):
            if not hgvs_g:
                continue
            snv = _snv_from_hgvs_g(hgvs_g, nc_to_chrom)
            if snv is None:
                continue
            hgvs_c = c_candidates[i] if i < len(c_candidates) else ""

            if revel_enabled and not (hgvs_c and revel_file_cache is not None and hgvs_c in revel_file_cache):
                revel_needed.add((
                    revel_mode, snv[0], snv[1], snv[2], snv[3],
                    revel_transcript_id or "", revel_aa_ref or "", revel_aa_alt or "",
                ))
            if am_enabled and not (hgvs_c and am_file_cache is not None and hgvs_c in am_file_cache):
                am_needed.add(snv)
            if mutpred2_dbnsfp_enabled:
                mp2_needed.add(snv)

    return revel_needed, am_needed, mp2_needed


def _prefetch_into_cache(
    keys: set[tuple[str, int, str, str]],
    lookup_one: Any,
    max_workers: int,
) -> dict[tuple[str, int, str, str], Any]:
    """Run *lookup_one* over *keys* concurrently and return a key -> result dict."""
    if not keys:
        return {}
    keys_list = list(keys)
    results: dict[tuple[str, int, str, str], Any] = {}
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for key, value in zip(keys_list, executor.map(lookup_one, keys_list)):
            results[key] = value
    return results


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=(
            "Annotate rows with REVEL and/or AlphaMissense missense pathogenicity scores "
            "using tabix-indexed pre-computed files."
        )
    )
    p.add_argument("input_file", help="Input CSV/TSV file")
    p.add_argument("output_file", help="Output CSV/TSV file")
    p.add_argument(
        "--revel-file",
        default=os.environ.get("REVEL_FILE"),
        metavar="PATH",
        help=(
            "bgzipped, tabix-indexed REVEL TSV (revel_hg38.tsv.gz). "
            "See module docstring for preparation instructions. "
            "Defaults to REVEL_FILE env var."
        ),
    )
    p.add_argument(
        "--revel-mode",
        choices=["coordinate", "transcript", "aa"],
        default="transcript",
        help=(
            "How REVEL's one-row-per-overlapping-transcript layout is disambiguated "
            "at a given genomic position (default: transcript). 'coordinate': max "
            "score across every matching row regardless of transcript (legacy "
            "behavior; works with the plain 5-column revel_hg38.tsv.gz). "
            "'transcript': require the row's Ensembl transcript ID to match this "
            "variant's own transcript (mapped from --mapped-hgvs-c-transcript-col "
            "via --revel-mane-file or --revel-transcript-mapping-file). 'aa': "
            "require the row's aaref/aaalt to match this variant's own amino-acid "
            "ref/alt (from --mapped-hgvs-p-ref-col/--mapped-hgvs-p-alt-col) — a "
            "variant with aa_ref == aa_alt (synonymous) can never match, since "
            "REVEL's rows are always missense. 'transcript' and 'aa' both require "
            "the extended REVEL file layout — see module docstring."
        ),
    )
    p.add_argument(
        "--revel-mane-file",
        default=os.environ.get("REVEL_MANE_FILE"),
        metavar="PATH",
        help=(
            "NCBI MANE summary file (MANE.GRCh38.*.summary.txt[.gz]) used to map this "
            "variant's RefSeq transcript to Ensembl for --revel-mode transcript. "
            "Mutually exclusive with --revel-transcript-mapping-file. "
            "Defaults to REVEL_MANE_FILE env var."
        ),
    )
    p.add_argument(
        "--revel-transcript-mapping-file",
        default=os.environ.get("REVEL_TRANSCRIPT_MAPPING_FILE"),
        metavar="PATH",
        help=(
            "Custom two-column TSV/CSV (headers source_id, target_id) mapping RefSeq "
            "transcript accessions to Ensembl, used by --revel-mode transcript instead "
            "of --revel-mane-file. Defaults to REVEL_TRANSCRIPT_MAPPING_FILE env var."
        ),
    )
    p.add_argument(
        "--mapped-hgvs-c-transcript-col",
        default="mapped_hgvs_c_transcript",
        help=(
            "Input column with the bare RefSeq transcript accession (e.g. "
            "'NM_058216.3'), used by --revel-mode transcript "
            "(default: mapped_hgvs_c_transcript)"
        ),
    )
    p.add_argument(
        "--alphamissense-file",
        default=os.environ.get("ALPHAMISSENSE_FILE"),
        metavar="PATH",
        help=(
            "bgzipped, tabix-indexed AlphaMissense TSV (AlphaMissense_hg38.tsv.gz). "
            "Defaults to ALPHAMISSENSE_FILE env var."
        ),
    )
    p.add_argument(
        "--dbnsfp-file",
        default=os.environ.get("DBNSFP_FILE"),
        metavar="PATH",
        help=(
            "bgzipped, tabix-indexed dbNSFP GRCh38 variant file "
            "(e.g. dbNSFP5.3.1a_grch38.gz). "
            "Used to annotate MutPred2_score. "
            "Ignored if --mutpred2-properties-file is also given. "
            "Defaults to DBNSFP_FILE env var."
        ),
    )
    p.add_argument(
        "--mutpred2-properties-file",
        default=os.environ.get("MUTPRED2_PROPERTIES_FILE"),
        metavar="PATH",
        help=(
            "MaveDB MP2-properties CSV (optionally gzipped), e.g. "
            "data_frame_missense_variants_MP2_properties.csv.gz. Preferred source for "
            "mutpred2.score. Join columns depend on --mutpred2-properties-join-key: "
            "'genomic' (default) requires mavedb_variant_urn, Chrom, hg38_start, "
            "hg38_end, ref_allele, alt_allele, 'MutPred2 score'; 'gene-aa' requires "
            "gene_symbol, AA, 'MutPred2 score'. "
            "Takes precedence over --dbnsfp-file when both are given. "
            "Defaults to MUTPRED2_PROPERTIES_FILE env var."
        ),
    )
    p.add_argument(
        "--mutpred2-properties-join-key",
        choices=["genomic", "gene-aa"],
        default="genomic",
        help=(
            "How --mutpred2-properties-file is joined to the input. 'genomic' "
            "(default): per reverse-translation candidate via --variant-urn-col and "
            "the --mapped-hgvs-g-*-col columns, so mutpred2.score is pipe-aligned like "
            "revel.score. 'gene-aa': gene symbol + amino-acid substitution (e.g. "
            "'T2A') via --gene-symbol-col and the --mapped-hgvs-p-ref/-start/-alt-col "
            "columns — MutPred2 is protein-level, so this matches regardless of which "
            "DNA reverse-translation candidate is used; the score is computed once per "
            "row and then duplicated across every pipe slot, so mutpred2.score is still "
            "pipe-aligned like the genomic join."
        ),
    )
    p.add_argument(
        "--mutpred2-gene-aa-long-indels",
        choices=["annotate", "ignore"],
        default="annotate",
        help=(
            "With --mutpred2-properties-join-key gene-aa, whether DNA reverse-"
            "translation candidates whose genomic ref or alt allele (from "
            "--mapped-hgvs-g-ref-col/--mapped-hgvs-g-alt-col) exceeds 3bp still "
            "receive the row's mutpred2.score. 'annotate' (default): every candidate "
            "gets the score regardless of allele length. 'ignore': candidates with "
            "max(len(ref), len(alt)) > 3bp get an empty slot instead. No effect with "
            "--mutpred2-properties-join-key genomic or --dbnsfp-file."
        ),
    )
    p.add_argument(
        "--mutpred2-gene-symbol-map-file",
        default=os.environ.get("MUTPRED2_GENE_SYMBOL_MAP_FILE"),
        metavar="PATH",
        help=(
            "With --mutpred2-properties-join-key gene-aa: a two-column TSV "
            "(gene_symbol, mutpred2_gene_symbol) remapping input gene symbols "
            "(from --gene-symbol-col) to the symbol used in "
            "--mutpred2-properties-file, for cases where a gene is named "
            "differently between the two (e.g. a combined-target dataset's "
            "'CALM1_2_3' vs. 'CALM1'). Only genes that differ need an entry; "
            "genes absent from the map are looked up unchanged. No effect with "
            "--mutpred2-properties-join-key genomic. Defaults to "
            "MUTPRED2_GENE_SYMBOL_MAP_FILE env var."
        ),
    )
    p.add_argument(
        "--variant-urn-col",
        default="variant_urn",
        help=(
            "Input column with the MaveDB variant URN, used as part of the lookup key "
            "for --mutpred2-properties-file with --mutpred2-properties-join-key genomic "
            "(default: variant_urn)"
        ),
    )
    p.add_argument(
        "--mapped-hgvs-g-chromosome-col",
        default="mapped_hgvs_g_chromosome",
        help=(
            "Input column with pipe-delimited genomic chromosome(s), used as part of "
            "the lookup key for --mutpred2-properties-file and --revel-training-file "
            "(default: mapped_hgvs_g_chromosome)"
        ),
    )
    p.add_argument(
        "--mapped-hgvs-g-start-col",
        default="mapped_hgvs_g_start",
        help=(
            "Input column with pipe-delimited genomic start position(s), used as part "
            "of the lookup key for --mutpred2-properties-file and --revel-training-file "
            "(default: mapped_hgvs_g_start)"
        ),
    )
    p.add_argument(
        "--mapped-hgvs-g-stop-col",
        default="mapped_hgvs_g_stop",
        help=(
            "Input column with pipe-delimited genomic stop position(s), used as part "
            "of the lookup key for --mutpred2-properties-file and --revel-training-file "
            "(default: mapped_hgvs_g_stop)"
        ),
    )
    p.add_argument(
        "--mapped-hgvs-g-ref-col",
        default="mapped_hgvs_g_ref",
        help=(
            "Input column with pipe-delimited genomic ref allele(s), used as part of "
            "the lookup key for --mutpred2-properties-file and --revel-training-file "
            "(default: mapped_hgvs_g_ref)"
        ),
    )
    p.add_argument(
        "--mapped-hgvs-g-alt-col",
        default="mapped_hgvs_g_alt",
        help=(
            "Input column with pipe-delimited genomic alt allele(s), used as part of "
            "the lookup key for --mutpred2-properties-file and --revel-training-file "
            "(default: mapped_hgvs_g_alt)"
        ),
    )
    p.add_argument(
        "--revel-training-file",
        default=os.environ.get("REVEL_TRAINING_FILE"),
        metavar="PATH",
        help=(
            "TSV of REVEL training variants with columns chromosome, hg38_start, "
            "hg38_end, ref_allele, alt_allele (e.g. data/revel_training_variants.tsv). "
            "Produces revel.train, looked up per reverse-translation candidate via the "
            "--mapped-hgvs-g-*-col columns, pipe-aligned like revel.score. "
            "Defaults to REVEL_TRAINING_FILE env var."
        ),
    )
    p.add_argument(
        "--mutpred2-training-file",
        default=os.environ.get("MUTPRED2_TRAINING_FILE"),
        metavar="PATH",
        help=(
            "TSV of MutPred2 training variants (e.g. data/mutpred2_training_variants.tsv), "
            "either with an 'hgvs_p' column (transcript/protein-accession-qualified, "
            "gene symbols ignored) or with 'gene_symbol' and 'unqualified_hgvs_p' columns. "
            "Produces mutpred2.train, evaluated once per row via --gene-symbol-col and "
            "--mapped-hgvs-p-col and duplicated across all DNA candidates. "
            "Defaults to MUTPRED2_TRAINING_FILE env var."
        ),
    )
    p.add_argument(
        "--gene-symbol-col",
        default="gene_symbol",
        help=(
            "Input column with the gene symbol, used as part of the lookup key for "
            "--mutpred2-training-file when that file has no 'hgvs_p' column, and for "
            "--mutpred2-properties-file with --mutpred2-properties-join-key gene-aa "
            "(default: gene_symbol)"
        ),
    )
    p.add_argument(
        "--mapped-hgvs-p-col",
        default="mapped_hgvs_p",
        help=(
            "Input column with protein HGVS value(s), used as part of the lookup key "
            "for --mutpred2-training-file (default: mapped_hgvs_p)"
        ),
    )
    p.add_argument(
        "--mapped-hgvs-p-ref-col",
        default="mapped_hgvs_p_ref",
        help=(
            "Input column with the single-letter reference amino acid, used with "
            "--mapped-hgvs-p-start-col and --mapped-hgvs-p-alt-col to build the AA "
            "substitution key (e.g. 'T2A') for --mutpred2-properties-file with "
            "--mutpred2-properties-join-key gene-aa (default: mapped_hgvs_p_ref). "
            "Also used (with --mapped-hgvs-p-alt-col) for --revel-mode aa."
        ),
    )
    p.add_argument(
        "--mapped-hgvs-p-start-col",
        default="mapped_hgvs_p_start",
        help=(
            "Input column with the amino acid position, used as part of the AA "
            "substitution key for --mutpred2-properties-file with "
            "--mutpred2-properties-join-key gene-aa (default: mapped_hgvs_p_start)"
        ),
    )
    p.add_argument(
        "--mapped-hgvs-p-alt-col",
        default="mapped_hgvs_p_alt",
        help=(
            "Input column with the single-letter alternate amino acid, used as part "
            "of the AA substitution key for --mutpred2-properties-file with "
            "--mutpred2-properties-join-key gene-aa (default: mapped_hgvs_p_alt). "
            "Also used (with --mapped-hgvs-p-ref-col) for --revel-mode aa."
        ),
    )
    p.add_argument(
        "--revel-cache-file",
        default=os.environ.get("REVEL_CACHE_FILE"),
        metavar="PATH",
        help=(
            "Path to a two-column TSV (hgvs, revel.score) used as a file-based "
            "REVEL cache. Looked up before tabix. Defaults to REVEL_CACHE_FILE env var."
        ),
    )
    p.add_argument(
        "--alphamissense-cache-file",
        default=os.environ.get("ALPHAMISSENSE_CACHE_FILE"),
        metavar="PATH",
        help=(
            "Path to a three-column TSV (hgvs, alphamissense.pathogenicity, "
            "alphamissense.class) used as a file-based AlphaMissense cache. "
            "Looked up before tabix. Defaults to ALPHAMISSENSE_CACHE_FILE env var."
        ),
    )
    p.add_argument(
        "--mapped-hgvs-g-col",
        default="mapped_hgvs_g",
        help="Input column containing pipe-delimited genomic HGVS values (default: mapped_hgvs_g)",
    )
    p.add_argument(
        "--mapped-hgvs-c-col",
        default="mapped_hgvs_c",
        help=(
            "Input column containing pipe-delimited transcript HGVS values used as "
            "keys for file-based caches (default: mapped_hgvs_c)"
        ),
    )
    p.add_argument(
        "--skip",
        type=int,
        default=0,
        help="Number of data rows to skip before annotation (default: 0)",
    )
    p.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Maximum number of data rows to annotate",
    )
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO)",
    )
    p.add_argument(
        "--csv-field-size-limit",
        type=int,
        default=csv.field_size_limit(),
        metavar="BYTES",
        help="Maximum per-field character length for CSV/TSV parsing (default: %(default)s)",
    )
    p.add_argument(
        "--max-workers",
        type=int,
        default=1,
        metavar="N",
        help=(
            "Number of parallel tabix lookup threads for REVEL/AlphaMissense/dbNSFP "
            "(default: 1). Values >1 make a first pass over the input to collect all "
            "unique SNVs, look them up concurrently, then annotate in a second pass. "
            "No effect on --mutpred2-properties-file, --revel-training-file, or "
            "--mutpred2-training-file, which don't use tabix."
        ),
    )
    return p.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    args = _parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    csv.field_size_limit(args.csv_field_size_limit)

    revel_path = Path(args.revel_file) if args.revel_file else None
    am_path = Path(args.alphamissense_file) if args.alphamissense_file else None
    dbnsfp_path = Path(args.dbnsfp_file) if args.dbnsfp_file else None
    mutpred2_properties_file = args.mutpred2_properties_file

    revel_file_cache: Optional[dict[str, str]] = None
    if args.revel_cache_file:
        revel_file_cache = _load_revel_file_cache(args.revel_cache_file)

    am_file_cache: Optional[dict[str, tuple[str, str]]] = None
    if args.alphamissense_cache_file:
        am_file_cache = _load_alphamissense_file_cache(args.alphamissense_cache_file)

    revel_enabled = revel_path is not None or revel_file_cache is not None
    am_enabled = am_path is not None or am_file_cache is not None

    if (
        not revel_enabled
        and not am_enabled
        and dbnsfp_path is None
        and not mutpred2_properties_file
        and not args.revel_training_file
        and not args.mutpred2_training_file
    ):
        logger.error(
            "At least one of --revel-file, --revel-cache-file, --alphamissense-file, "
            "--alphamissense-cache-file, --mutpred2-properties-file, --dbnsfp-file, "
            "--revel-training-file, or --mutpred2-training-file must be provided."
        )
        raise SystemExit(1)

    if revel_path is not None and not revel_path.exists():
        logger.error("REVEL file not found: %s", revel_path)
        raise SystemExit(1)

    revel_transcript_mapping: dict[str, str] = {}
    revel_transcript_base_mapping: dict[str, str] = {}
    if revel_path is not None and args.revel_mode != "coordinate":
        try:
            _validate_revel_file_for_mode(revel_path, args.revel_mode)
        except ValueError as exc:
            logger.error(str(exc))
            raise SystemExit(1)
        if args.revel_mode == "transcript":
            if bool(args.revel_mane_file) == bool(args.revel_transcript_mapping_file):
                logger.error(
                    "--revel-mode transcript requires exactly one of --revel-mane-file "
                    "or --revel-transcript-mapping-file."
                )
                raise SystemExit(1)
            if args.revel_mane_file:
                if not Path(args.revel_mane_file).exists():
                    logger.error("REVEL MANE mapping file not found: %s", args.revel_mane_file)
                    raise SystemExit(1)
                revel_transcript_mapping = _build_mane_mapping(args.revel_mane_file, "to-ensembl")
            else:
                if not Path(args.revel_transcript_mapping_file).exists():
                    logger.error(
                        "REVEL transcript mapping file not found: %s",
                        args.revel_transcript_mapping_file,
                    )
                    raise SystemExit(1)
                revel_transcript_mapping = _build_custom_mapping(args.revel_transcript_mapping_file)
            revel_transcript_base_mapping = _build_base_mapping(revel_transcript_mapping)
    if am_path is not None and not am_path.exists():
        logger.error("AlphaMissense file not found: %s", am_path)
        raise SystemExit(1)
    if mutpred2_properties_file and not Path(mutpred2_properties_file).exists():
        logger.error("MutPred2 properties file not found: %s", mutpred2_properties_file)
        raise SystemExit(1)
    if args.revel_training_file and not Path(args.revel_training_file).exists():
        logger.error("REVEL training-variants file not found: %s", args.revel_training_file)
        raise SystemExit(1)
    if args.mutpred2_training_file and not Path(args.mutpred2_training_file).exists():
        logger.error("MutPred2 training-variants file not found: %s", args.mutpred2_training_file)
        raise SystemExit(1)
    if args.mutpred2_gene_symbol_map_file and not Path(args.mutpred2_gene_symbol_map_file).exists():
        logger.error(
            "MutPred2 gene-symbol-map file not found: %s", args.mutpred2_gene_symbol_map_file
        )
        raise SystemExit(1)

    mutpred2_properties_cache: Optional[dict[tuple[str, str, int, int, str, str], str]] = None
    mutpred2_gene_aa_cache: Optional[dict[tuple[str, str], str]] = None
    mutpred2_gene_symbol_map: Optional[dict[str, str]] = None
    if mutpred2_properties_file:
        if args.mutpred2_properties_join_key == "gene-aa":
            mutpred2_gene_aa_cache = _load_mutpred2_properties_gene_aa_cache(mutpred2_properties_file)
            if args.mutpred2_gene_symbol_map_file:
                mutpred2_gene_symbol_map = _load_mutpred2_gene_symbol_map(
                    args.mutpred2_gene_symbol_map_file
                )
        else:
            mutpred2_properties_cache = _load_mutpred2_properties_file_cache(mutpred2_properties_file)
        if dbnsfp_path is not None:
            logger.warning(
                "Both --dbnsfp-file and --mutpred2-properties-file given; using "
                "--mutpred2-properties-file for mutpred2.score and ignoring --dbnsfp-file."
            )
    if dbnsfp_path is not None and not dbnsfp_path.exists():
        logger.error("dbNSFP file not found: %s", dbnsfp_path)
        raise SystemExit(1)

    revel_training_set: Optional[set[tuple[str, int, int, str, str]]] = None
    if args.revel_training_file:
        revel_training_set = _load_revel_training_variants(args.revel_training_file)

    mutpred2_training_schema: Optional[str] = None
    mutpred2_training_set: Optional[set] = None
    if args.mutpred2_training_file:
        mutpred2_training_schema, mutpred2_training_set = _load_mutpred2_training_variants(
            args.mutpred2_training_file
        )

    # Only require tabix when at least one tabix-indexed file is configured.
    if revel_path is not None or am_path is not None or dbnsfp_path is not None:
        if subprocess.run(["tabix", "--version"], capture_output=True, check=False).returncode not in (0, 1):
            logger.error("tabix executable not found; install htslib.")
            raise SystemExit(1)

    input_path = Path(args.input_file)
    output_path = Path(args.output_file)
    delim = "\t" if input_path.suffix.lower() in (".tsv", ".txt") else ","

    mutpred2_enabled = (
        mutpred2_properties_cache is not None
        or mutpred2_gene_aa_cache is not None
        or dbnsfp_path is not None
    )

    ann_cols: list[str] = []
    if revel_enabled:
        ann_cols.extend(REVEL_COLS)
    if am_enabled:
        ann_cols.extend(ALPHAMISSENSE_COLS)
    if mutpred2_enabled:
        ann_cols.extend(DBNSFP_COLS)
    if revel_training_set is not None:
        ann_cols.extend(REVEL_TRAIN_COLS)
    if mutpred2_training_set is not None:
        ann_cols.extend(MUTPRED2_TRAIN_COLS)

    revel_cache: dict[tuple, Optional[str]] = {}
    am_cache: dict[tuple[str, int, str, str], Optional[tuple[str, str]]] = {}
    mutpred2_cache: dict[tuple[str, int, str, str], Optional[str]] = {}

    mutpred2_dbnsfp_enabled = (
        dbnsfp_path is not None
        and mutpred2_properties_cache is None
        and mutpred2_gene_aa_cache is None
    )
    if args.max_workers > 1 and (revel_path is not None or am_path is not None or mutpred2_dbnsfp_enabled):
        logger.info("Prefetch pass: scanning input for unique SNVs (max_workers=%d)…", args.max_workers)
        with input_path.open("r", encoding="utf-8", newline="") as scan_fh:
            scan_reader = csv.DictReader(scan_fh, delimiter=delim)
            scan_rows = islice(
                scan_reader,
                args.skip,
                None if args.limit is None else args.skip + args.limit,
            )
            revel_needed, am_needed, mp2_needed = _collect_needed_snv_keys(
                scan_rows,
                nc_to_chrom=NC_TO_CHROM_GRCH38,
                mapped_hgvs_g_col=args.mapped_hgvs_g_col,
                mapped_hgvs_c_col=args.mapped_hgvs_c_col,
                revel_enabled=revel_path is not None,
                revel_file_cache=revel_file_cache,
                am_enabled=am_path is not None,
                am_file_cache=am_file_cache,
                mutpred2_dbnsfp_enabled=mutpred2_dbnsfp_enabled,
                revel_mode=args.revel_mode,
                revel_transcript_mapping=revel_transcript_mapping,
                revel_transcript_base_mapping=revel_transcript_base_mapping,
                mapped_hgvs_c_transcript_col=args.mapped_hgvs_c_transcript_col,
                mapped_hgvs_p_ref_col=args.mapped_hgvs_p_ref_col,
                mapped_hgvs_p_alt_col=args.mapped_hgvs_p_alt_col,
            )

        logger.info(
            "Prefetch pass: %d unique REVEL / %d AlphaMissense / %d MutPred2(dbNSFP) SNVs to look up",
            len(revel_needed), len(am_needed), len(mp2_needed),
        )
        prefetch_started = time.monotonic()
        if revel_needed:
            revel_cache.update(_prefetch_into_cache(
                revel_needed,
                lambda key: _lookup_revel(
                    revel_path, key[1], key[2], key[3], key[4], {},  # type: ignore[arg-type]
                    mode=key[0],
                    transcript_id=key[5] or None,
                    aa_ref=key[6] or None,
                    aa_alt=key[7] or None,
                ),
                args.max_workers,
            ))
        if am_needed:
            am_cache.update(_prefetch_into_cache(
                am_needed,
                lambda key: _lookup_alphamissense(am_path, key[0], key[1], key[2], key[3], {}),  # type: ignore[arg-type]
                args.max_workers,
            ))
        if mp2_needed:
            mutpred2_cache.update(_prefetch_into_cache(
                mp2_needed,
                lambda key: _lookup_mutpred2(dbnsfp_path, key[0], key[1], key[2], key[3], {}),  # type: ignore[arg-type]
                args.max_workers,
            ))
        logger.info("Prefetch pass complete in %.1fs.", time.monotonic() - prefetch_started)

    processed = 0
    scored_revel = 0
    scored_am = 0
    scored_mutpred2 = 0
    revel_train_rows_matched = 0
    revel_train_variants_matched = 0
    revel_train_variants_total = 0
    mutpred2_train_rows_matched = 0
    mutpred2_train_variants_matched = 0
    mutpred2_train_variants_total = 0
    started = time.monotonic()

    with input_path.open("r", encoding="utf-8", newline="") as in_fh, \
         output_path.open("w", encoding="utf-8", newline="") as out_fh:

        reader = csv.DictReader(in_fh, delimiter=delim)
        if reader.fieldnames is None:
            logger.error("Input file appears empty: %s", input_path)
            raise SystemExit(1)

        fieldnames = list(reader.fieldnames)
        out_fieldnames = fieldnames + [c for c in ann_cols if c not in fieldnames]
        writer = csv.DictWriter(
            out_fh,
            fieldnames=out_fieldnames,
            delimiter=delim,
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()

        selected = islice(
            reader,
            args.skip,
            None if args.limit is None else args.skip + args.limit,
        )

        for row in selected:
            ann = annotate_row(
                row,
                nc_to_chrom=NC_TO_CHROM_GRCH38,
                mapped_hgvs_g_col=args.mapped_hgvs_g_col,
                mapped_hgvs_c_col=args.mapped_hgvs_c_col,
                revel_path=revel_path,
                alphamissense_path=am_path,
                dbnsfp_path=dbnsfp_path,
                revel_cache=revel_cache,
                revel_mode=args.revel_mode,
                revel_transcript_mapping=revel_transcript_mapping,
                revel_transcript_base_mapping=revel_transcript_base_mapping,
                mapped_hgvs_c_transcript_col=args.mapped_hgvs_c_transcript_col,
                am_cache=am_cache,
                mutpred2_cache=mutpred2_cache,
                revel_file_cache=revel_file_cache,
                am_file_cache=am_file_cache,
                mutpred2_properties_cache=mutpred2_properties_cache,
                mutpred2_gene_aa_cache=mutpred2_gene_aa_cache,
                mutpred2_gene_aa_long_indels=args.mutpred2_gene_aa_long_indels,
                mutpred2_gene_symbol_map=mutpred2_gene_symbol_map,
                variant_urn_col=args.variant_urn_col,
                mapped_hgvs_g_chromosome_col=args.mapped_hgvs_g_chromosome_col,
                mapped_hgvs_g_start_col=args.mapped_hgvs_g_start_col,
                mapped_hgvs_g_stop_col=args.mapped_hgvs_g_stop_col,
                mapped_hgvs_g_ref_col=args.mapped_hgvs_g_ref_col,
                mapped_hgvs_g_alt_col=args.mapped_hgvs_g_alt_col,
                revel_training_set=revel_training_set,
                mutpred2_training_schema=mutpred2_training_schema,
                mutpred2_training_set=mutpred2_training_set,
                gene_symbol_col=args.gene_symbol_col,
                mapped_hgvs_p_col=args.mapped_hgvs_p_col,
                mapped_hgvs_p_ref_col=args.mapped_hgvs_p_ref_col,
                mapped_hgvs_p_start_col=args.mapped_hgvs_p_start_col,
                mapped_hgvs_p_alt_col=args.mapped_hgvs_p_alt_col,
            )
            row.update(ann)
            writer.writerow(row)

            processed += 1
            if revel_enabled and row.get("revel.score"):
                scored_revel += 1
            if am_enabled and row.get("alphamissense.pathogenicity"):
                scored_am += 1
            if mutpred2_enabled and row.get("mutpred2.score"):
                scored_mutpred2 += 1
            if revel_training_set is not None:
                vals = ann["revel.train"].split("|")
                n_true = vals.count("true")
                revel_train_variants_total += len(vals)
                revel_train_variants_matched += n_true
                if n_true:
                    revel_train_rows_matched += 1
            if mutpred2_training_set is not None:
                vals = ann["mutpred2.train"].split("|")
                n_true = vals.count("true")
                mutpred2_train_variants_total += len(vals)
                mutpred2_train_variants_matched += n_true
                if n_true:
                    mutpred2_train_rows_matched += 1

            if processed % 1000 == 0:
                elapsed = max(time.monotonic() - started, 1e-9)
                logger.info(
                    "Progress: %d rows processed (%.1f rows/s)",
                    processed,
                    processed / elapsed,
                )

    elapsed = max(time.monotonic() - started, 1e-9)
    if revel_path is not None:
        logger.info(
            "REVEL: %d/%d rows scored (mode: %s, cache: %d unique SNVs queried)",
            scored_revel, processed, args.revel_mode, len(revel_cache),
        )
    if am_path is not None:
        logger.info(
            "AlphaMissense: %d/%d rows scored (cache: %d unique SNVs queried)",
            scored_am, processed, len(am_cache),
        )
    if mutpred2_properties_cache is not None:
        logger.info(
            "MutPred2 (properties file, genomic join): %d/%d rows scored "
            "(%d candidates in file cache)",
            scored_mutpred2, processed, len(mutpred2_properties_cache),
        )
    elif mutpred2_gene_aa_cache is not None:
        logger.info(
            "MutPred2 (properties file, gene-AA join): %d/%d rows scored "
            "(%d (gene, AA) keys in file cache); joined on %s/%s+%s+%s vs. "
            "properties-file gene_symbol/AA; long indels (>3bp ref/alt): %s; "
            "gene-symbol map: %s",
            scored_mutpred2, processed, len(mutpred2_gene_aa_cache),
            args.gene_symbol_col, args.mapped_hgvs_p_ref_col,
            args.mapped_hgvs_p_start_col, args.mapped_hgvs_p_alt_col,
            args.mutpred2_gene_aa_long_indels,
            f"{len(mutpred2_gene_symbol_map)} entries" if mutpred2_gene_symbol_map else "none",
        )
    elif dbnsfp_path is not None:
        logger.info(
            "MutPred2 (dbNSFP): %d/%d rows scored (cache: %d unique SNVs queried)",
            scored_mutpred2, processed, len(mutpred2_cache),
        )
    if revel_training_set is not None:
        logger.info(
            "REVEL training-set overlap: %d/%d assayed variants matched "
            "(%d/%d genomic variants); joined on %s/%s/%s/%s/%s vs. REVEL-training "
            "chromosome/hg38_start/hg38_end/ref_allele/alt_allele",
            revel_train_rows_matched, processed,
            revel_train_variants_matched, revel_train_variants_total,
            args.mapped_hgvs_g_chromosome_col, args.mapped_hgvs_g_start_col,
            args.mapped_hgvs_g_stop_col, args.mapped_hgvs_g_ref_col, args.mapped_hgvs_g_alt_col,
        )
    if mutpred2_training_set is not None:
        if mutpred2_training_schema == "qualified":
            join_desc = f"{args.mapped_hgvs_p_col} vs. MutPred2-training hgvs_p (qualified; gene symbols ignored)"
        else:
            join_desc = (
                f"{args.gene_symbol_col}/{args.mapped_hgvs_p_col} vs. "
                "MutPred2-training gene_symbol/unqualified_hgvs_p"
            )
        logger.info(
            "MutPred2 training-set overlap: %d/%d assayed variants matched "
            "(%d/%d genomic variants); joined on %s",
            mutpred2_train_rows_matched, processed,
            mutpred2_train_variants_matched, mutpred2_train_variants_total,
            join_desc,
        )
    logger.info(
        "Done. %d rows written to %s (%.1f rows/s)",
        processed, output_path, processed / elapsed,
    )


if __name__ == "__main__":
    main()
