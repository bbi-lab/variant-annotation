"""backfill_clingen_allele_id.py

Provenance
----------
This module is project-specific, but its ClinGen lookup helpers and column
selection behavior are intentionally aligned with src.map_variants, which in
turn mirrors MaveDB behavior. If src.map_variants is treated as AGPL-coupled,
review this file together with it before choosing a more permissive license.

Populate the ``clingen_allele_id`` column (or a custom-named column) for rows
where it is currently blank, by querying the ClinGen Allele Registry using the
already-mapped HGVS strings in the file.

This script backfills the *protein-level* ``clingen_allele_id`` column written
by Step 1 (``map_variants``).  It is not a replacement for Step 3
(``add_dna_clingen_allele_ids``), which resolves *DNA-level* allele IDs for
reverse-translated candidates.  Use this script when ``map_variants`` left some
``clingen_allele_id`` cells blank due to transient API errors or rate-limiting
during the initial run.

HGVS priority (matches map_variants.py behaviour):
  1. mapped_hgvs_c  (c. / n. notation)
  2. mapped_hgvs_g  (g. notation)
  3. mapped_hgvs_p  (p. notation)

The script selects the first non-empty value from the three HGVS columns and
submits it as a single HGVS string to the ClinGen Allele Registry.  If the
column contains a pipe-delimited list (i.e. the file has already been through
Step 2 reverse translation), the entire pipe-delimited string would be used
as the query, which will fail.  Apply this script only to Step-1 output.

Blank-node-style response identifiers (``_:PA...`` / ``_:CA...``) that
ClinGen emits when no real allele record exists are treated as misses and
left blank.

Rows that already have a value in the clingen_allele_id column are left
untouched.  Row order in the output is preserved.

Usage
-----
    python src/backfill_clingen_allele_id.py input.tsv output.tsv [options]

Options
-------
    --clingen-allele-id COLNAME   Column name for the allele ID
                                  (default: clingen_allele_id)
    --hgvs-c-col COLNAME          Column containing mapped c./n. HGVS
                                  (default: mapped_hgvs_c)
    --hgvs-g-col COLNAME          Column containing mapped g. HGVS
                                  (default: mapped_hgvs_g)
    --hgvs-p-col COLNAME          Column containing mapped p. HGVS
                                  (default: mapped_hgvs_p)
    --concurrency N               Max concurrent ClinGen requests (default: 5)
    --max-retries N               Retries per request (default: 3)
    --skip N                      Skip the first N data rows
    --limit N                     Process at most N rows (after --skip)
    --log-level LEVEL             Logging verbosity (default: INFO)
"""

from __future__ import annotations

import argparse
import asyncio
import csv
import logging
import sys
import time
from pathlib import Path
from typing import Optional

import requests
from dotenv import load_dotenv

load_dotenv()

logger = logging.getLogger(__name__)

CLINGEN_API_URL = "https://reg.genome.network/allele"
CLINGEN_RETRY_DELAY = 2.0
PROGRESS_EVERY_ROWS = 1000


# ---------------------------------------------------------------------------
# ClinGen helpers (adapted from map_variants.py)
# ---------------------------------------------------------------------------


def _detect_separator(file_path: str) -> str:
    return "\t" if Path(file_path).suffix.lower() in (".tsv", ".txt") else ","


def _query_clingen_by_hgvs(hgvs_string: str, max_retries: int = 3) -> Optional[dict]:
    for attempt in range(max_retries):
        try:
            resp = requests.get(
                CLINGEN_API_URL,
                params={"hgvs": hgvs_string},
                timeout=30,
                headers={"Accept": "application/json"},
            )
            if resp.status_code == 200:
                return resp.json()
            if resp.status_code == 404:
                logger.warning("ClinGen 404 for %s", hgvs_string)
                return None
            if resp.status_code == 429:
                wait = CLINGEN_RETRY_DELAY * (2**attempt)
                logger.warning(
                    "ClinGen rate-limited for %s; waiting %.1f s.", hgvs_string, wait
                )
                time.sleep(wait)
                continue
            logger.warning(
                "ClinGen returned HTTP %d for %s", resp.status_code, hgvs_string
            )
            return None
        except requests.exceptions.RequestException as exc:
            logger.warning(
                "ClinGen request failed for %s (attempt %d/%d): %s",
                hgvs_string,
                attempt + 1,
                max_retries,
                exc,
            )
            if attempt < max_retries - 1:
                time.sleep(CLINGEN_RETRY_DELAY)
    return None


async def _query_clingen_batch(
    hgvs_strings: list[str],
    max_concurrency: int,
    max_retries: int,
) -> dict[str, Optional[dict]]:
    semaphore = asyncio.Semaphore(max_concurrency)

    async def _query_one(hgvs: str) -> tuple[str, Optional[dict]]:
        async with semaphore:
            loop = asyncio.get_event_loop()
            result = await loop.run_in_executor(
                None, _query_clingen_by_hgvs, hgvs, max_retries
            )
            return hgvs, result

    tasks = [_query_one(h) for h in hgvs_strings]
    pairs = await asyncio.gather(*tasks)
    return {hgvs: data for hgvs, data in pairs}


def _extract_clingen_allele_id(data: dict) -> Optional[str]:
    """Extract ClinGen allele identifier (for example ``CA123456``) from response.

    Returns None for placeholder values such as ``_:PA...`` or ``_:CA...`` that
    ClinGen emits when no real allele record exists.
    """
    def _is_real_id(value: str) -> bool:
        # Reject blank-node-style identifiers like "_:PA..." or "_:CA..."
        return bool(value) and not value.startswith("_:")

    at_id: str = data.get("@id", "") or ""
    if at_id:
        fragment = at_id.rstrip("/").rsplit("/", 1)[-1]
        return fragment if _is_real_id(fragment) else None

    fallback = data.get("id")
    if isinstance(fallback, str):
        stripped = fallback.strip()
        if _is_real_id(stripped):
            return stripped
    return None


# ---------------------------------------------------------------------------
# Main backfill logic
# ---------------------------------------------------------------------------


def _pick_hgvs(
    row: dict,
    hgvs_c_col: str,
    hgvs_g_col: str,
    hgvs_p_col: str,
) -> Optional[str]:
    """Return the best available HGVS string from *row* (c > g > p)."""
    for col in (hgvs_c_col, hgvs_g_col, hgvs_p_col):
        val = (row.get(col) or "").strip()
        if val:
            return val
    return None


def backfill(
    input_path: str,
    output_path: str,
    *,
    clingen_allele_id_col: str = "clingen_allele_id",
    hgvs_c_col: str = "mapped_hgvs_c",
    hgvs_g_col: str = "mapped_hgvs_g",
    hgvs_p_col: str = "mapped_hgvs_p",
    max_concurrency: int = 5,
    max_retries: int = 3,
    skip: int = 0,
    limit: Optional[int] = None,
) -> None:
    sep = _detect_separator(input_path)
    out_sep = _detect_separator(output_path)

    with open(input_path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter=sep)
        if reader.fieldnames is None:
            raise ValueError(f"Input file {input_path!r} appears to be empty.")
        fieldnames: list[str] = list(reader.fieldnames)
        rows: list[dict] = list(reader)

    if skip:
        rows = rows[skip:]
    if limit is not None:
        rows = rows[:limit]

    # Ensure the output column exists
    if clingen_allele_id_col not in fieldnames:
        fieldnames.append(clingen_allele_id_col)
        logger.info(
            "Column %r not found in input; it will be appended.", clingen_allele_id_col
        )

    # Identify rows that need a lookup
    needs_lookup: list[int] = []
    hgvs_for_row: dict[int, str] = {}

    for idx, row in enumerate(rows):
        existing = (row.get(clingen_allele_id_col) or "").strip()
        if existing:
            continue
        hgvs = _pick_hgvs(row, hgvs_c_col, hgvs_g_col, hgvs_p_col)
        if hgvs:
            needs_lookup.append(idx)
            hgvs_for_row[idx] = hgvs

    logger.info(
        "%d of %d rows need ClinGen lookup.",
        len(needs_lookup),
        len(rows),
    )

    if not needs_lookup:
        logger.info("Nothing to do; writing output unchanged.")
    else:
        # Deduplicate HGVS strings to minimise API calls
        unique_hgvs: list[str] = list(dict.fromkeys(hgvs_for_row[i] for i in needs_lookup))
        logger.info(
            "Querying ClinGen for %d unique HGVS string(s) (concurrency=%d)…",
            len(unique_hgvs),
            max_concurrency,
        )

        results: dict[str, Optional[dict]] = asyncio.run(
            _query_clingen_batch(unique_hgvs, max_concurrency, max_retries)
        )

        resolved = 0
        for idx in needs_lookup:
            hgvs = hgvs_for_row[idx]
            data = results.get(hgvs)
            allele_id: Optional[str] = None
            if data is not None:
                allele_id = _extract_clingen_allele_id(data)
            rows[idx][clingen_allele_id_col] = allele_id or ""
            if allele_id:
                resolved += 1
            if idx < 10 or idx % PROGRESS_EVERY_ROWS == 0:
                logger.debug(
                    "Row %d: hgvs=%r allele_id=%r", idx, hgvs, allele_id
                )

        logger.info(
            "Resolved allele ID for %d of %d rows that needed lookup.",
            resolved,
            len(needs_lookup),
        )

    with open(output_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter=out_sep,
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

    logger.info("Wrote %d rows to %s", len(rows), output_path)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Backfill clingen_allele_id for rows where it is blank.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input", help="Input CSV/TSV file path.")
    p.add_argument("output", help="Output CSV/TSV file path.")
    p.add_argument(
        "--clingen-allele-id",
        dest="clingen_allele_id_col",
        default="clingen_allele_id",
        metavar="COLNAME",
        help="Column name for the ClinGen allele ID.",
    )
    p.add_argument(
        "--hgvs-c-col",
        default="mapped_hgvs_c",
        metavar="COLNAME",
        help="Column containing mapped c./n. HGVS strings.",
    )
    p.add_argument(
        "--hgvs-g-col",
        default="mapped_hgvs_g",
        metavar="COLNAME",
        help="Column containing mapped g. HGVS strings.",
    )
    p.add_argument(
        "--hgvs-p-col",
        default="mapped_hgvs_p",
        metavar="COLNAME",
        help="Column containing mapped p. HGVS strings.",
    )
    p.add_argument(
        "--concurrency",
        type=int,
        default=5,
        metavar="N",
        help="Max concurrent ClinGen API requests.",
    )
    p.add_argument(
        "--max-retries",
        type=int,
        default=3,
        metavar="N",
        help="Retries per ClinGen request.",
    )
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity level.",
    )
    p.add_argument(
        "--csv-field-size-limit",
        type=int,
        default=csv.field_size_limit(),
        metavar="BYTES",
        help="Maximum per-field character length for CSV/TSV parsing (default: %(default)s).",
    )
    p.add_argument(
        "--skip",
        type=int,
        default=0,
        metavar="N",
        help="Number of data rows to skip from the start of the input.",
    )
    p.add_argument(
        "--limit",
        type=int,
        default=None,
        metavar="N",
        help="Maximum number of data rows to process after applying --skip.",
    )
    return p


def main(argv: Optional[list[str]] = None) -> None:
    args = _build_parser().parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )
    csv.field_size_limit(args.csv_field_size_limit)

    backfill(
        args.input,
        args.output,
        clingen_allele_id_col=args.clingen_allele_id_col,
        hgvs_c_col=args.hgvs_c_col,
        hgvs_g_col=args.hgvs_g_col,
        hgvs_p_col=args.hgvs_p_col,
        max_concurrency=args.concurrency,
        max_retries=args.max_retries,
        skip=args.skip,
        limit=args.limit,
    )


if __name__ == "__main__":
    main()
