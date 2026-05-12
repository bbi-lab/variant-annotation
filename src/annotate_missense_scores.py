"""Annotate variants with pre-computed missense pathogenicity scores.

Supported scores (hg38 / GRCh38 only):

  REVEL (Rare Exome Variant Ensemble Learner)
  --revel-file revel_hg38.tsv.gz
    Source: https://sites.google.com/site/revelgenomics/downloads
    Range:  0–1 (higher = more likely pathogenic)
    Scope:  missense SNVs only

  AlphaMissense (Google DeepMind)
  --alphamissense-file AlphaMissense_hg38.tsv.gz
    Source: https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
    Index:  generated locally with ``tabix -s 1 -b 2 -e 2 AlphaMissense_hg38.tsv.gz``
    Range:  0–1 (higher = more likely pathogenic)
    Scope:  missense SNVs only

  MutPred2 (via dbNSFP)
  --dbnsfp-file dbNSFP5.3.1a_grch38.gz
    Source: https://sites.google.com/site/jpopgen/dbNSFP
    The file and its .tbi index are available pre-built (see run script for URLs).
    Range:  0–1 (higher = more likely pathogenic)
    Scope:  missense SNVs only; returns max score across transcripts

At least one of --revel-file, --alphamissense-file, or --dbnsfp-file is required.


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
                   {print $1"\\t"$3"\\t"$4"\\t"$5"\\t"$8}' \\
      | (printf '#chr\\tpos\\tref\\talt\\trevel_score\\n'; sort -k1,1V -k2,2n) \\
      | bgzip > revel_hg38.tsv.gz
    tabix -s 1 -b 2 -e 2 -S 1 revel_hg38.tsv.gz

The resulting file has five tab-separated columns::

    #chr  pos  ref  alt  revel_score


Output columns
--------------

  revel.score                — REVEL score string (empty for non-SNV / non-missense)
  alphamissense.pathogenicity — AlphaMissense pathogenicity score (0–1)
  alphamissense.class        — likely_benign / ambiguous / likely_pathogenic
  mutpred2.score             — MutPred2 score from dbNSFP (single value per row; protein-level model)

For rows with pipe-delimited genomic HGVS values the REVEL and AlphaMissense
columns are pipe-aligned to match the input candidate positions.  MutPred2 is a
protein-level model, so all reverse-translation candidates encode the same amino
acid substitution; a single score (the maximum across candidates) is emitted.
"""

from __future__ import annotations

import argparse
import csv
from itertools import islice
import logging
import os
import re
import subprocess
import time
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)

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


def _run_tabix(path: Path, chrom: str, pos: int) -> list[str]:
    """Return non-comment lines from a tabix point query at *chrom*:*pos*."""
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

def _lookup_revel(
    path: Path,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    cache: dict[tuple[str, int, str, str], Optional[str]],
) -> Optional[str]:
    """Return the maximum REVEL score string for an SNV, or ``None`` if absent.

    REVEL may have multiple rows per position (one per transcript); we return
    the maximum score across all matching rows.

    Expected prepared-file column layout (tab-separated, 0-indexed):
      0: chr  1: pos  2: ref  3: alt  4: revel_score
    """
    key = (chrom, pos, ref, alt)
    if key in cache:
        return cache[key]

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
# Row-level annotation
# ---------------------------------------------------------------------------

def annotate_row(
    row: dict[str, str],
    *,
    nc_to_chrom: dict[str, str],
    mapped_hgvs_g_col: str,
    revel_path: Optional[Path],
    alphamissense_path: Optional[Path],
    dbnsfp_path: Optional[Path] = None,
    revel_cache: dict[tuple[str, int, str, str], Optional[str]],
    am_cache: dict[tuple[str, int, str, str], Optional[tuple[str, str]]],
    mutpred2_cache: Optional[dict[tuple[str, int, str, str], Optional[str]]] = None,
) -> dict[str, str]:
    """Return annotation columns for a single row.

    Output values are pipe-aligned to the pipe-delimited candidates in the
    ``mapped_hgvs_g_col`` input column.  Non-SNV candidates produce empty
    strings in every score column.
    """
    candidates = _split_pipe((row.get(mapped_hgvs_g_col) or "").strip())

    revel_vals: list[str] = []
    am_path_vals: list[str] = []
    am_class_vals: list[str] = []
    mutpred2_vals: list[str] = []

    # Use a local throwaway cache if the caller didn't supply one (preserves
    # backward-compatibility; cross-row deduplication requires a real dict).
    _mp2_cache: dict[tuple[str, int, str, str], Optional[str]] = (
        mutpred2_cache if mutpred2_cache is not None else {}
    )

    for hgvs in candidates:
        snv = _snv_from_hgvs_g(hgvs, nc_to_chrom) if hgvs else None

        if snv is not None and revel_path is not None:
            chrom, pos, ref, alt = snv
            r = _lookup_revel(revel_path, chrom, pos, ref, alt, revel_cache)
            revel_vals.append(r or "")
        else:
            revel_vals.append("")

        if snv is not None and alphamissense_path is not None:
            chrom, pos, ref, alt = snv
            a = _lookup_alphamissense(alphamissense_path, chrom, pos, ref, alt, am_cache)
            if a is not None:
                am_path_vals.append(a[0])
                am_class_vals.append(a[1])
            else:
                am_path_vals.append("")
                am_class_vals.append("")
        else:
            am_path_vals.append("")
            am_class_vals.append("")

        if snv is not None and dbnsfp_path is not None:
            chrom, pos, ref, alt = snv
            m = _lookup_mutpred2(dbnsfp_path, chrom, pos, ref, alt, _mp2_cache)
            if m is not None:
                mutpred2_vals.append(m)
        # (non-SNV candidates are simply skipped for the protein-level score)

    sep = "|" if len(candidates) > 1 else ""
    out: dict[str, str] = {}
    if revel_path is not None:
        out["revel.score"] = sep.join(revel_vals)
    if alphamissense_path is not None:
        out["alphamissense.pathogenicity"] = sep.join(am_path_vals)
        out["alphamissense.class"] = sep.join(am_class_vals)
    if dbnsfp_path is not None:
        # MutPred2 is protein-level: all candidates encode the same amino acid
        # substitution, so emit the single best score across candidates.
        if mutpred2_vals:
            best_mp2 = max(mutpred2_vals, key=float)
        else:
            best_mp2 = ""
        out["mutpred2.score"] = best_mp2
    return out


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
            "Defaults to DBNSFP_FILE env var."
        ),
    )
    p.add_argument(
        "--mapped-hgvs-g-col",
        default="mapped_hgvs_g",
        help="Input column containing pipe-delimited genomic HGVS values (default: mapped_hgvs_g)",
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

    if revel_path is None and am_path is None and dbnsfp_path is None:
        logger.error(
            "At least one of --revel-file, --alphamissense-file, or --dbnsfp-file "
            "must be provided."
        )
        raise SystemExit(1)

    if revel_path is not None and not revel_path.exists():
        logger.error("REVEL file not found: %s", revel_path)
        raise SystemExit(1)
    if am_path is not None and not am_path.exists():
        logger.error("AlphaMissense file not found: %s", am_path)
        raise SystemExit(1)
    if dbnsfp_path is not None and not dbnsfp_path.exists():
        logger.error("dbNSFP file not found: %s", dbnsfp_path)
        raise SystemExit(1)

    # Check that tabix is available.
    if subprocess.run(["tabix", "--version"], capture_output=True, check=False).returncode not in (0, 1):
        logger.error("tabix executable not found; install htslib.")
        raise SystemExit(1)

    input_path = Path(args.input_file)
    output_path = Path(args.output_file)
    delim = "\t" if input_path.suffix.lower() in (".tsv", ".txt") else ","

    ann_cols: list[str] = []
    if revel_path is not None:
        ann_cols.extend(REVEL_COLS)
    if am_path is not None:
        ann_cols.extend(ALPHAMISSENSE_COLS)
    if dbnsfp_path is not None:
        ann_cols.extend(DBNSFP_COLS)

    revel_cache: dict[tuple[str, int, str, str], Optional[str]] = {}
    am_cache: dict[tuple[str, int, str, str], Optional[tuple[str, str]]] = {}
    mutpred2_cache: dict[tuple[str, int, str, str], Optional[str]] = {}

    processed = 0
    scored_revel = 0
    scored_am = 0
    scored_mutpred2 = 0
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
                revel_path=revel_path,
                alphamissense_path=am_path,
                dbnsfp_path=dbnsfp_path,
                revel_cache=revel_cache,
                am_cache=am_cache,
                mutpred2_cache=mutpred2_cache,
            )
            row.update(ann)
            writer.writerow(row)

            processed += 1
            if revel_path is not None and row.get("revel.score"):
                scored_revel += 1
            if am_path is not None and row.get("alphamissense.pathogenicity"):
                scored_am += 1
            if dbnsfp_path is not None and row.get("mutpred2.score"):
                scored_mutpred2 += 1

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
            "REVEL: %d/%d rows scored (cache: %d unique SNVs queried)",
            scored_revel, processed, len(revel_cache),
        )
    if am_path is not None:
        logger.info(
            "AlphaMissense: %d/%d rows scored (cache: %d unique SNVs queried)",
            scored_am, processed, len(am_cache),
        )
    if dbnsfp_path is not None:
        logger.info(
            "MutPred2 (dbNSFP): %d/%d rows scored (cache: %d unique SNVs queried)",
            scored_mutpred2, processed, len(mutpred2_cache),
        )
    logger.info(
        "Done. %d rows written to %s (%.1f rows/s)",
        processed, output_path, processed / elapsed,
    )


if __name__ == "__main__":
    main()
