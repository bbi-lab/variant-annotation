"""Add a DNA-level ClinGen allele-ID column.

This script writes one output column containing DNA-level ClinGen allele IDs:

- DNA rows: a single ID equal to ``clingen_allele_id`` (when present).
- Protein reverse-translation rows: one ID per reverse-translated DNA candidate,
  aligned to the pipe-delimited ``mapped_hgvs_c`` / ``mapped_hgvs_g`` values.

For each DNA candidate pair, ClinGen is queried in this order:
  1. HGVS c. string
  2. HGVS g. string (fallback if c. is absent or not found)

If a candidate cannot be resolved, an empty value is kept in that position. The
resulting column is pipe-delimited and preserves candidate cardinality, e.g.
``CA1||CA3``.

Usage:
    python -m src.add_dna_clingen_allele_ids input.tsv output.tsv
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import logging
from threading import Lock
from pathlib import Path
from typing import Optional

from src.lib.clingen import query_clingen_by_hgvs

logger = logging.getLogger(__name__)

def _detect_separator(file_path: str) -> str:
    return "\t" if Path(file_path).suffix.lower() in (".tsv", ".txt") else ","


def _query_clingen_by_hgvs(hgvs_string: str, max_retries: int = 3) -> Optional[dict]:
    """Query ClinGen Allele Registry by HGVS string."""
    return query_clingen_by_hgvs(hgvs_string, max_retries=max_retries)


def _extract_clingen_allele_id(data: dict) -> Optional[str]:
    def _normalize_clingen_id(value: str) -> str:
        """Return empty string for placeholder ClinGen IDs (e.g. '_:CA...', '_:PA...')."""
        text = (value or "").strip()
        if not text:
            return ""
        if text.startswith("_:"):
            return ""
        return text

    at_id: str = data.get("@id", "") or ""
    if at_id:
        fragment = at_id.rstrip("/").rsplit("/", 1)[-1]
        normalized = _normalize_clingen_id(fragment)
        return normalized or None
    fallback = data.get("id")
    if isinstance(fallback, str) and fallback.strip():
        normalized = _normalize_clingen_id(fallback)
        return normalized or None
    return None


def _split_pipe(value: str) -> list[str]:
    """Split a pipe-delimited string while preserving empty entries."""
    if "|" not in value:
        return [value.strip()]
    return [part.strip() for part in value.split("|")]


def _candidate_pairs(hgvs_c: str, hgvs_g: str) -> list[tuple[str, str]]:
    """Return aligned DNA-candidate pairs from mapped_hgvs_c and mapped_hgvs_g.

    If neither string contains a pipe, this yields one pair.
    If either contains a pipe, values are treated as candidate lists and aligned
    by index. Length mismatches are padded with empty strings.
    """
    c_val = (hgvs_c or "").strip()
    g_val = (hgvs_g or "").strip()

    if not c_val and not g_val:
        return []

    has_multi = "|" in c_val or "|" in g_val
    if not has_multi:
        return [(c_val, g_val)]

    c_parts = _split_pipe(c_val) if c_val else [""]
    g_parts = _split_pipe(g_val) if g_val else [""]

    if len(c_parts) != len(g_parts):
        logger.warning(
            "Mismatched candidate counts: mapped_hgvs_c=%d mapped_hgvs_g=%d; padding shorter side",
            len(c_parts),
            len(g_parts),
        )

    n = max(len(c_parts), len(g_parts))
    c_parts.extend([""] * (n - len(c_parts)))
    g_parts.extend([""] * (n - len(g_parts)))
    return list(zip(c_parts, g_parts))


def _lookup_allele_id_for_candidate(
    hgvs_c: str,
    hgvs_g: str,
    *,
    max_retries: int,
    lookup_cache: dict[str, str],
    cache_lock: Optional[Lock] = None,
) -> str:
    """Resolve one DNA candidate to a ClinGen allele ID using c-then-g lookup."""
    for hgvs in (hgvs_c, hgvs_g):
        query = (hgvs or "").strip()
        if not query:
            continue

        if cache_lock is None:
            if query not in lookup_cache:
                data = _query_clingen_by_hgvs(query, max_retries=max_retries)
                lookup_cache[query] = (_extract_clingen_allele_id(data) or "") if data else ""
            allele_id = lookup_cache.get(query, "")
        else:
            with cache_lock:
                cached = lookup_cache.get(query)
            if cached is None:
                data = _query_clingen_by_hgvs(query, max_retries=max_retries)
                resolved = (_extract_clingen_allele_id(data) or "") if data else ""
                with cache_lock:
                    # Keep first resolved value if another thread wrote it first.
                    lookup_cache.setdefault(query, resolved)
                    cached = lookup_cache.get(query, "")
            allele_id = cached or ""

        if allele_id:
            return allele_id
    return ""


def _is_dna_variant_row(
    row: dict[str, str],
    *,
    raw_hgvs_nt_col: str,
    raw_hgvs_pro_col: str,
) -> bool:
    """Return True when a row is clearly DNA-based (not protein-only)."""
    raw_nt = (row.get(raw_hgvs_nt_col) or "").strip()
    # raw_pro = (row.get(raw_hgvs_pro_col) or "").strip()
    return bool(raw_nt) # and not bool(raw_pro)


def _validate_clingen_id_prefix(
    existing_id: str,
    is_dna_row: bool,
    row_index: int,
) -> None:
    """Validate that ClinGen allele ID has the correct prefix based on variant type.
    
    DNA variants should have IDs starting with "CA" (ClinGen Allele).
    Protein variants should have IDs starting with "PA" (Protein Allele).
    
    Raises ValueError if the prefix doesn't match the variant type.
    """
    if not existing_id:
        return
    
    expected_prefix = "CA" if is_dna_row else "PA"
    if not existing_id.startswith(expected_prefix):
        variant_type = "DNA" if is_dna_row else "protein"
        raise ValueError(
            f"Row {row_index}: Invalid ClinGen allele ID prefix for {variant_type} variant. "
            f"Expected prefix '{expected_prefix}' but got '{existing_id[:2]}' in '{existing_id}'"
        )


def build_dna_clingen_value(
    row: dict[str, str],
    *,
    clingen_allele_id_col: str,
    hgvs_c_col: str,
    hgvs_g_col: str,
    raw_hgvs_nt_col: str,
    raw_hgvs_pro_col: str,
    max_retries: int,
    lookup_cache: dict[str, str],
    cache_lock: Optional[Lock] = None,
) -> str:
    """Build the DNA-level ClinGen allele-ID cell value for one row."""
    existing_id = (row.get(clingen_allele_id_col) or "").strip()
    if existing_id.startswith("_:"):
        existing_id = ""
    pairs = _candidate_pairs(row.get(hgvs_c_col, ""), row.get(hgvs_g_col, ""))

    if not pairs:
        return ""

    # Reuse existing clingen_allele_id only for explicit DNA rows.
    if (
        len(pairs) == 1
        and existing_id
        and "|" not in existing_id
        and _is_dna_variant_row(
            row,
            raw_hgvs_nt_col=raw_hgvs_nt_col,
            raw_hgvs_pro_col=raw_hgvs_pro_col,
        )
    ):
        return existing_id

    ids = [
        _lookup_allele_id_for_candidate(
            c,
            g,
            max_retries=max_retries,
            lookup_cache=lookup_cache,
            cache_lock=cache_lock,
        )
        for c, g in pairs
    ]
    return "|".join(ids)


def _process_row_for_dna_ids(
    row_idx: int,
    row: dict[str, str],
    *,
    clingen_allele_id_col: str,
    hgvs_c_col: str,
    hgvs_g_col: str,
    raw_hgvs_nt_col: str,
    raw_hgvs_pro_col: str,
    max_retries: int,
    output_col: str,
    lookup_cache: dict[str, str],
    cache_lock: Lock,
) -> tuple[int, dict[str, str], bool]:
    """Compute DNA ClinGen value for one row and return ordered result tuple."""
    row_index_for_errors = row_idx + 2  # Account for header row.
    existing_id = (row.get(clingen_allele_id_col) or "").strip()
    if existing_id.startswith("_:"):
        existing_id = ""
    is_dna_row = _is_dna_variant_row(
        row,
        raw_hgvs_nt_col=raw_hgvs_nt_col,
        raw_hgvs_pro_col=raw_hgvs_pro_col,
    )

    _validate_clingen_id_prefix(existing_id, is_dna_row, row_index_for_errors)

    value = build_dna_clingen_value(
        row,
        clingen_allele_id_col=clingen_allele_id_col,
        hgvs_c_col=hgvs_c_col,
        hgvs_g_col=hgvs_g_col,
        raw_hgvs_nt_col=raw_hgvs_nt_col,
        raw_hgvs_pro_col=raw_hgvs_pro_col,
        max_retries=max_retries,
        lookup_cache=lookup_cache,
        cache_lock=cache_lock,
    )
    row[output_col] = value
    return row_idx, row, bool(value)


def add_dna_clingen_allele_ids(
    input_path: str,
    output_path: str,
    *,
    output_col: str = "dna_clingen_allele_id",
    clingen_allele_id_col: str = "clingen_allele_id",
    hgvs_c_col: str = "mapped_hgvs_c",
    hgvs_g_col: str = "mapped_hgvs_g",
    raw_hgvs_nt_col: str = "raw_hgvs_nt",
    raw_hgvs_pro_col: str = "raw_hgvs_pro",
    max_retries: int = 3,
    max_workers: int = 8,
) -> None:
    """Read input table, add DNA-level ClinGen ID column, write output table.

    Rows are written in input order as soon as contiguous results are available,
    while ClinGen lookups are performed concurrently across rows.

    Rows are submitted concurrently and written in input order when ready.
    """
    if max_workers < 1:
        raise ValueError("max_workers must be >= 1")

    in_sep = _detect_separator(input_path)
    out_sep = _detect_separator(output_path)

    lookup_cache: dict[str, str] = {}
    cache_lock = Lock()
    populated = 0
    written = 0

    with open(input_path, newline="", encoding="utf-8") as in_fh, open(
        output_path, "w", newline="", encoding="utf-8"
    ) as out_fh:
        reader = csv.DictReader(in_fh, delimiter=in_sep)
        if reader.fieldnames is None:
            raise ValueError(f"Input file {input_path!r} appears to be empty.")
        fieldnames: list[str] = list(reader.fieldnames)

        if output_col not in fieldnames:
            fieldnames.append(output_col)

        writer = csv.DictWriter(
            out_fh,
            fieldnames=fieldnames,
            delimiter=out_sep,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        out_fh.flush()

        pending_by_index: dict[int, tuple[dict[str, str], bool]] = {}
        next_to_write = 0

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {
                executor.submit(
                    _process_row_for_dna_ids,
                    idx,
                    row,
                    clingen_allele_id_col=clingen_allele_id_col,
                    hgvs_c_col=hgvs_c_col,
                    hgvs_g_col=hgvs_g_col,
                    raw_hgvs_nt_col=raw_hgvs_nt_col,
                    raw_hgvs_pro_col=raw_hgvs_pro_col,
                    max_retries=max_retries,
                    output_col=output_col,
                    lookup_cache=lookup_cache,
                    cache_lock=cache_lock,
                ): idx
                for idx, row in enumerate(reader)
            }

            for fut in as_completed(futures):
                row_idx, row_out, has_value = fut.result()
                pending_by_index[row_idx] = (row_out, has_value)

                # Stream contiguous rows immediately while preserving order.
                while next_to_write in pending_by_index:
                    ready_row, ready_has_value = pending_by_index.pop(next_to_write)
                    writer.writerow(ready_row)
                    out_fh.flush()
                    written += 1
                    if ready_has_value:
                        populated += 1
                    next_to_write += 1

    logger.info(
        "Wrote %d rows to %s (%d rows with non-empty %s)",
        written,
        output_path,
        populated,
        output_col,
    )


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=(
            "Add a DNA-level ClinGen allele ID column. Protein reverse-translation "
            "rows produce a pipe-delimited ID list aligned with mapped_hgvs_c/g candidates."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("input", help="Input CSV/TSV file path")
    p.add_argument("output", help="Output CSV/TSV file path")
    p.add_argument(
        "--output-col",
        default="dna_clingen_allele_id",
        help="Output column name for DNA-level ClinGen allele IDs",
    )
    p.add_argument(
        "--clingen-allele-id-col",
        default="clingen_allele_id",
        help="Column containing existing ClinGen allele ID from map_variants",
    )
    p.add_argument(
        "--hgvs-c-col",
        default="mapped_hgvs_c",
        help="Column containing mapped HGVS c strings",
    )
    p.add_argument(
        "--hgvs-g-col",
        default="mapped_hgvs_g",
        help="Column containing mapped HGVS g strings",
    )
    p.add_argument(
        "--raw-hgvs-nt-col",
        default="raw_hgvs_nt",
        help="Column containing original DNA HGVS input",
    )
    p.add_argument(
        "--raw-hgvs-pro-col",
        default="raw_hgvs_pro",
        help="Column containing original protein HGVS input",
    )
    p.add_argument(
        "--max-retries",
        type=int,
        default=3,
        help="Retries per ClinGen request",
    )
    p.add_argument(
        "--max-workers",
        type=int,
        default=8,
        help="Concurrent worker threads for ClinGen lookups",
    )
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity",
    )
    return p


def main(argv: Optional[list[str]] = None) -> None:
    args = _build_parser().parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
        datefmt="%H:%M:%S",
    )

    add_dna_clingen_allele_ids(
        args.input,
        args.output,
        output_col=args.output_col,
        clingen_allele_id_col=args.clingen_allele_id_col,
        hgvs_c_col=args.hgvs_c_col,
        hgvs_g_col=args.hgvs_g_col,
        raw_hgvs_nt_col=args.raw_hgvs_nt_col,
        raw_hgvs_pro_col=args.raw_hgvs_pro_col,
        max_retries=args.max_retries,
        max_workers=args.max_workers,
    )


if __name__ == "__main__":
    main()
