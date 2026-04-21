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
is tried in order and the first successful ClinVar lookup is used.

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
import csv
import gzip
import io
import logging
import os
import sys
import time
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Optional

import requests

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

CLINGEN_API_URL = "https://reg.genome.network/allele"

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

# Retry settings for the ClinGen API
CLINGEN_MAX_RETRIES = 3
CLINGEN_RETRY_DELAY = 2.0  # seconds between retries

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
# ClinGen → ClinVar allele ID resolution
# ---------------------------------------------------------------------------


def resolve_clinvar_allele_id(
    clingen_id: str,
    cache: dict[str, str],
    max_retries: int = CLINGEN_MAX_RETRIES,
    retry_delay: float = CLINGEN_RETRY_DELAY,
) -> str:
    """Return the ClinVar allele ID for *clingen_id*, or empty string if not found.

    Results are stored in *cache* (in-memory throughout the run) to avoid redundant
    API calls when the same ClinGen ID appears multiple times.

    Args:
        clingen_id: A ClinGen Allele Registry ID, e.g. ``CA123456``.
        cache: Mutable dict used as an in-process cache (modified in-place).
        max_retries: Number of retry attempts on transient HTTP errors.
        retry_delay: Seconds to wait between retries.

    Returns:
        ClinVar allele ID string (integer as string, e.g. ``"12345"``), or ``""`` if
        no association exists or the allele is not registered in ClinGen.
    """
    if clingen_id in cache:
        return cache[clingen_id]

    url = f"{CLINGEN_API_URL}/{urllib.parse.quote(clingen_id, safe='')}"
    last_exc: Optional[Exception] = None

    for attempt in range(max_retries):
        try:
            response = requests.get(url, timeout=30)
            if response.status_code == 404:
                cache[clingen_id] = ""
                return ""
            response.raise_for_status()
            api_data = response.json()
            allele_id = (
                api_data.get("externalRecords", {})
                .get("ClinVarAlleles", [{}])[0]
                .get("alleleId")
            )
            result = str(allele_id) if allele_id is not None else ""
            cache[clingen_id] = result
            return result
        except (requests.HTTPError, requests.ConnectionError, requests.Timeout) as exc:
            last_exc = exc
            if isinstance(exc, requests.HTTPError) and exc.response is not None:
                status = exc.response.status_code
                # Don't retry client errors other than rate-limiting
                if status != 429 and 400 <= status < 500:
                    logger.warning(
                        "ClinGen API returned %s for %s; skipping", status, clingen_id
                    )
                    cache[clingen_id] = ""
                    return ""
            if attempt < max_retries - 1:
                logger.warning(
                    "ClinGen API error for %s (attempt %d/%d): %s; retrying in %.1fs",
                    clingen_id,
                    attempt + 1,
                    max_retries,
                    exc,
                    retry_delay,
                )
                time.sleep(retry_delay)

    logger.error(
        "ClinGen API failed for %s after %d attempts: %s",
        clingen_id,
        max_retries,
        last_exc,
    )
    cache[clingen_id] = ""
    return ""


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
    out = {
        f"{col_prefix}.clinical_significance": "",
        f"{col_prefix}.review_status": "",
        f"{col_prefix}.stars": "",
        f"{col_prefix}.last_review_date": "",
    }

    raw_clingen = row.get(dna_clingen_allele_id_col, "").strip()
    if not raw_clingen:
        return out

    # Support pipe-delimited candidate IDs (from reverse_translate_protein_variants)
    candidates = [c.strip() for c in raw_clingen.split("|") if c.strip()]

    for clingen_id in candidates:
        clinvar_id = resolve_clinvar_allele_id(clingen_id, clingen_cache)
        if not clinvar_id:
            continue
        record = clinvar_data.get(clinvar_id)
        if record is None:
            continue
        sig = record.get(COL_CLINICAL_SIGNIFICANCE, "")
        review = record.get(COL_REVIEW_STATUS, "")
        last_eval = record.get(COL_LAST_EVALUATED, "")
        out[f"{col_prefix}.clinical_significance"] = sig
        out[f"{col_prefix}.review_status"] = review
        out[f"{col_prefix}.stars"] = stars_for_review_status(review)
        out[f"{col_prefix}.last_review_date"] = last_eval
        return out

    return out


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
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    args = _parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )

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
        for row in reader:
            annotations = annotate_row(
                row,
                clinvar_data,
                clingen_cache,
                col_prefix,
                dna_clingen_allele_id_col=args.dna_clingen_allele_id_col,
            )
            row.update(annotations)
            writer.writerow(row)
            processed += 1
            if annotations[f"{col_prefix}.clinical_significance"]:
                annotated += 1
            if processed % 500 == 0:
                logger.info("Processed %d rows (%d annotated) …", processed, annotated)

    logger.info(
        "Done. %d rows processed, %d annotated with ClinVar data → %s",
        processed,
        annotated,
        out_path,
    )


if __name__ == "__main__":
    main()
