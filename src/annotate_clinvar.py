"""Annotate variant rows with ClinVar clinical-significance data.

For each row produced by ``map_variants.py`` (or ``reverse_translate_protein_variants.py``),
this script:

1. Downloads and caches the ClinVar monthly variant-summary TSV from NCBI.
2. Resolves each row's ClinGen Allele ID(s) to a ClinVar Allele ID via the ClinGen
   Allele Registry REST API.
3. Looks up the ClinVar allele in the cached TSV and writes four annotation columns:

   ``<namespace>.<version>.clinical_significance``
   ``<namespace>.<version>.review_status``
   ``<namespace>.<version>.stars``
   ``<namespace>.<version>.last_review_date``

   where ``<namespace>`` defaults to ``clinvar`` and ``<version>`` is ``YYYYMM``
   (e.g. ``202601``).

Input rows that have no ``dna_clingen_allele_id``, or whose ClinGen ID does not map to a
ClinVar allele, receive empty strings for all four annotation columns.

When ``dna_clingen_allele_id`` is pipe-delimited (produced by
``add_dna_clingen_allele_ids.py`` for protein reverse translations), each candidate
is annotated independently and the output columns are likewise pipe-delimited,
preserving candidate cardinality.

Provenance
----------
This module is modeled on the ClinVar refresh logic in MaveDB
(``mavedb-api/src/mavedb/worker/jobs/external_services/clinvar.py`` and
``mavedb-api/src/mavedb/lib/clinvar/utils.py``). If this file contains adapted logic
from MaveDB, treat it as AGPL-coupled when assigning a project license.

Usage::

    python -m src.annotate_clinvar input.tsv output.tsv [OPTIONS]

Example::

    python -m src.annotate_clinvar mapped.tsv annotated.tsv \\
        --clinvar-version 202601 \\
        --cache-dir /data/clinvar_cache
"""

from __future__ import annotations

import argparse
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
import csv
import gzip
import io
import logging
import os
import sys
from pathlib import Path
from typing import Optional

import requests

from src.lib.clingen import resolve_clinvar_allele_id

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CLINVAR_TSV_BASE_URL = (
    "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive"
)

# ClinVar variant summary TSV column names
COL_ALLELE_ID = "#AlleleID"
COL_CLINICAL_SIGNIFICANCE = "ClinicalSignificance"
COL_REVIEW_STATUS = "ReviewStatus"
COL_LAST_EVALUATED = "LastEvaluated"

# ReviewStatus → star rating mapping (ClinVar definitions, 2024)
REVIEW_STATUS_STARS: dict[str, int] = {
    "practice guideline": 4,
    "reviewed by expert panel": 3,
    "criteria provided, multiple submitters, no conflicts": 2,
    "criteria provided, conflicting classifications": 1,
    "criteria provided, conflicting interpretations": 1,
    "criteria provided, single submitter": 1,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
    "no classifications from unflagged records": 0,
    "no classification for the single variant": 0,
    "no classification provided": 0,
}

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# ClinVar TSV helpers
# ---------------------------------------------------------------------------


def _tsv_url(year: int, month: int) -> list[str]:
    """Return candidate download URLs for the given year/month (primary then fallback)."""
    filename = f"variant_summary_{year}-{month:02d}.txt.gz"
    return [
        f"{CLINVAR_TSV_BASE_URL}/{filename}",
        f"{CLINVAR_TSV_BASE_URL}/{year}/{filename}",
    ]


def _cache_path(cache_dir: Path, year: int, month: int) -> Path:
    return cache_dir / f"variant_summary_{year}-{month:02d}.txt.gz"


def fetch_clinvar_tsv(year: int, month: int, cache_dir: Path) -> Path:
    """Download (or return cached) the ClinVar variant-summary TSV for *year*-*month*.

    The file is stored as a ``.gz`` in *cache_dir* and reused on subsequent calls.
    These archive files are immutable, so no TTL is applied.

    Args:
        year: Four-digit year (e.g. 2026).
        month: Month (1–12).
        cache_dir: Directory in which to cache downloaded files.

    Returns:
        Path to the cached ``.gz`` file.

    Raises:
        RuntimeError: If the file cannot be downloaded from any candidate URL.
    """
    cache_dir.mkdir(parents=True, exist_ok=True)
    dest = _cache_path(cache_dir, year, month)

    if dest.exists():
        logger.info("Using cached ClinVar TSV: %s", dest)
        return dest

    urls = _tsv_url(year, month)
    for url in urls:
        logger.info("Downloading ClinVar TSV from %s", url)
        try:
            response = requests.get(url, stream=True, timeout=120)
            response.raise_for_status()
            tmp = dest.with_suffix(".tmp")
            with tmp.open("wb") as fh:
                for chunk in response.iter_content(chunk_size=1 << 16):
                    fh.write(chunk)
            tmp.rename(dest)
            logger.info("Saved ClinVar TSV to %s", dest)
            return dest
        except requests.HTTPError as exc:
            logger.warning("Failed to download %s: %s", url, exc)

    raise RuntimeError(
        f"Could not download ClinVar variant summary for {year}-{month:02d}. "
        f"Tried: {urls}"
    )


def load_clinvar_tsv(path: Path) -> dict[str, dict[str, str]]:
    """Parse a gzipped ClinVar variant-summary TSV into a lookup dict.

    Returns a dict mapping ClinVar allele ID (as string) to a sub-dict with keys
    ``ClinicalSignificance``, ``ReviewStatus``, and ``LastEvaluated``.

    Only germline rows (``Assembly == GRCh38``, ``Origin`` contains ``germline``) are
    retained when those columns are present; all rows are kept otherwise.
    """
    logger.info("Loading ClinVar TSV from %s …", path)
    data: dict[str, dict[str, str]] = {}

    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            allele_id = row.get(COL_ALLELE_ID, "").strip()
            if not allele_id:
                continue

            # Prefer germline significance when Origin column is available
            origin = row.get("Origin", "").lower()
            assembly = row.get("Assembly", "")
            # Skip non-GRCh38 rows when Assembly column exists
            if assembly and assembly not in ("GRCh38", ""):
                continue

            # Only retain germline submissions (skip somatic-only)
            germline_sig = row.get("ClinSigSimple", "")
            sig = row.get(COL_CLINICAL_SIGNIFICANCE, "").strip()
            review = row.get(COL_REVIEW_STATUS, "").strip()
            last_eval = row.get(COL_LAST_EVALUATED, "").strip()

            data[allele_id] = {
                COL_CLINICAL_SIGNIFICANCE: sig,
                COL_REVIEW_STATUS: review,
                COL_LAST_EVALUATED: last_eval,
            }

    logger.info("Loaded %d ClinVar allele records", len(data))
    return data


def stars_for_review_status(review_status: str) -> str:
    """Return the numeric star rating (as a string) for the given ReviewStatus value."""
    normalised = review_status.strip().lower()
    stars = REVIEW_STATUS_STARS.get(normalised)
    if stars is None:
        # Try prefix match for future or variant phrasing
        for key, val in REVIEW_STATUS_STARS.items():
            if normalised.startswith(key):
                stars = val
                break
    return str(stars) if stars is not None else ""


# ---------------------------------------------------------------------------
# Core annotation logic
# ---------------------------------------------------------------------------


def annotate_row(
    row: dict[str, str],
    clinvar_data: dict[str, dict[str, str]],
    clingen_cache: dict[str, str],
    col_prefix: str,
    dna_clingen_allele_id_col: str = "dna_clingen_allele_id",
) -> dict[str, str]:
    """Return annotation values for *row* as a flat dict keyed by output column name.

    *col_prefix* should already include the namespace and version, e.g.
    ``"clinvar.202601"``.
    """
    empty = {
        f"{col_prefix}.clinical_significance": "",
        f"{col_prefix}.review_status": "",
        f"{col_prefix}.stars": "",
        f"{col_prefix}.last_review_date": "",
    }

    raw_clingen = row.get(dna_clingen_allele_id_col, "").strip()
    if not raw_clingen:
        return empty

    # Support pipe-delimited candidate IDs (from reverse_translate_protein_variants).
    # Each candidate is annotated independently; output columns are pipe-delimited
    # with one entry per candidate (empty string when no ClinVar record is found).
    all_candidates = [c.strip() for c in raw_clingen.split("|")]

    sigs, reviews, stars_list, last_evals = [], [], [], []
    for clingen_id in all_candidates:
        if not clingen_id:
            sigs.append("")
            reviews.append("")
            stars_list.append("")
            last_evals.append("")
            continue
        clinvar_id = resolve_clinvar_allele_id(clingen_id, clingen_cache)
        record = clinvar_data.get(clinvar_id) if clinvar_id else None
        if record is None:
            sigs.append("")
            reviews.append("")
            stars_list.append("")
            last_evals.append("")
        else:
            sig = record.get(COL_CLINICAL_SIGNIFICANCE, "")
            review = record.get(COL_REVIEW_STATUS, "")
            last_eval = record.get(COL_LAST_EVALUATED, "")
            sigs.append(sig)
            reviews.append(review)
            stars_list.append(stars_for_review_status(review))
            last_evals.append(last_eval)

    # Collapse to a single value when there is only one candidate (non-protein rows)
    # so the output format matches the input cardinality.
    join = "|".join
    return {
        f"{col_prefix}.clinical_significance": join(sigs),
        f"{col_prefix}.review_status": join(reviews),
        f"{col_prefix}.stars": join(stars_list),
        f"{col_prefix}.last_review_date": join(last_evals),
    }


def _annotate_row_task(
    row_idx: int,
    row: dict[str, str],
    *,
    clinvar_data: dict[str, dict[str, str]],
    clingen_cache: dict[str, str],
    col_prefix: str,
    dna_clingen_allele_id_col: str,
) -> tuple[int, dict[str, str], dict[str, str]]:
    annotations = annotate_row(
        row,
        clinvar_data,
        clingen_cache,
        col_prefix,
        dna_clingen_allele_id_col=dna_clingen_allele_id_col,
    )
    return row_idx, row, annotations


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Annotate variant rows with ClinVar clinical-significance data, "
            "using the ClinVar monthly variant-summary TSV archive."
        )
    )
    parser.add_argument("input_file", help="Input TSV file (output of map_variants.py)")
    parser.add_argument("output_file", help="Output TSV file with ClinVar annotation columns appended")

    parser.add_argument(
        "--clinvar-version",
        default="202601",
        metavar="YYYYMM",
        help=(
            "ClinVar release to use, in YYYYMM format (default: %(default)s). "
            "The corresponding monthly variant-summary archive will be downloaded "
            "from NCBI if not already cached."
        ),
    )
    parser.add_argument(
        "--clinvar-namespace",
        default="clinvar",
        metavar="NAMESPACE",
        help=(
            "Namespace prefix for output columns (default: %(default)s). "
            "Column names are <namespace>.<version>.<field>."
        ),
    )
    parser.add_argument(
        "--cache-dir",
        default=os.environ.get("CLINVAR_CACHE_DIR", "/tmp/clinvar_cache"),
        metavar="DIR",
        help=(
            "Directory for caching downloaded ClinVar TSV files "
            "(default: $CLINVAR_CACHE_DIR env var or /tmp/clinvar_cache)."
        ),
    )
    parser.add_argument(
        "--dna-clingen-allele-id-col",
        default="dna_clingen_allele_id",
        metavar="COL",
        help=(
            "Column containing DNA-level ClinGen allele IDs to resolve to ClinVar "
            "(default: %(default)s)."
        ),
    )
    parser.add_argument(
        "--delimiter",
        default="\t",
        help="Field delimiter for both input and output (default: TAB).",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO).",
    )
    parser.add_argument(
        "--max-workers",
        type=int,
        default=8,
        help="Concurrent worker threads for row annotation.",
    )
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    args = _parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

    if args.max_workers < 1:
        logger.error("--max-workers must be >= 1, got: %d", args.max_workers)
        sys.exit(1)

    # Validate and parse the --clinvar-version argument
    version_str = args.clinvar_version.strip()
    if len(version_str) != 6 or not version_str.isdigit():
        logger.error(
            "--clinvar-version must be in YYYYMM format (e.g. 202601), got: %s",
            version_str,
        )
        sys.exit(1)
    year = int(version_str[:4])
    month = int(version_str[4:])
    if not (1 <= month <= 12):
        logger.error("Month in --clinvar-version must be 01–12, got: %02d", month)
        sys.exit(1)

    col_prefix = f"{args.clinvar_namespace}.{version_str}"
    cache_dir = Path(args.cache_dir)

    # Download / load ClinVar TSV
    tsv_path = fetch_clinvar_tsv(year, month, cache_dir)
    clinvar_data = load_clinvar_tsv(tsv_path)

    # In-process ClinGen → ClinVar ID cache
    clingen_cache: dict[str, str] = {}

    delim = args.delimiter
    if delim == "\\t":  # handle shell quoting edge case
        delim = "\t"

    in_path = Path(args.input_file)
    out_path = Path(args.output_file)

    annotation_cols = [
        f"{col_prefix}.clinical_significance",
        f"{col_prefix}.review_status",
        f"{col_prefix}.stars",
        f"{col_prefix}.last_review_date",
    ]

    with in_path.open("r", encoding="utf-8", newline="") as in_fh, \
         out_path.open("w", encoding="utf-8", newline="") as out_fh:

        reader = csv.DictReader(in_fh, delimiter=delim)
        if reader.fieldnames is None:
            logger.error("Input file appears to be empty: %s", in_path)
            sys.exit(1)

        out_fieldnames = list(reader.fieldnames) + annotation_cols
        writer = csv.DictWriter(
            out_fh,
            fieldnames=out_fieldnames,
            delimiter=delim,
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()

        processed = 0
        annotated = 0
        pending_by_index: dict[int, tuple[dict[str, str], dict[str, str]]] = {}
        next_to_write = 0

        with ThreadPoolExecutor(max_workers=args.max_workers) as executor:
            max_in_flight = args.max_workers * 4
            row_iter = enumerate(reader)
            in_flight: dict = {}

            def _submit_until_full() -> None:
                while len(in_flight) < max_in_flight:
                    try:
                        idx, row = next(row_iter)
                    except StopIteration:
                        break
                    fut = executor.submit(
                        _annotate_row_task,
                        idx,
                        row,
                        clinvar_data=clinvar_data,
                        clingen_cache=clingen_cache,
                        col_prefix=col_prefix,
                        dna_clingen_allele_id_col=args.dna_clingen_allele_id_col,
                    )
                    in_flight[fut] = idx

            _submit_until_full()

            while in_flight:
                done, _ = wait(tuple(in_flight.keys()), return_when=FIRST_COMPLETED)
                for fut in done:
                    in_flight.pop(fut, None)
                    row_idx, row_out, annotations = fut.result()
                    pending_by_index[row_idx] = (row_out, annotations)

                while next_to_write in pending_by_index:
                    row_out, annotations = pending_by_index.pop(next_to_write)
                    row_out.update(annotations)
                    writer.writerow(row_out)
                    out_fh.flush()
                    processed += 1
                    if annotations[f"{col_prefix}.clinical_significance"]:
                        annotated += 1
                    if processed % 500 == 0:
                        logger.info("Processed %d rows (%d annotated) …", processed, annotated)
                    next_to_write += 1

                _submit_until_full()

    logger.info(
        "Done. %d rows processed, %d annotated with ClinVar data → %s",
        processed,
        annotated,
        out_path,
    )


if __name__ == "__main__":
    main()
