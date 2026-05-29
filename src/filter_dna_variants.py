#!/usr/bin/env python3
"""Filter DNA candidate lists in place while preserving row count.

This script keeps one input row per output row, but removes DNA candidates that
are not in the requested retained classes. The default behavior keeps SNVs and
wild-type/unchanged candidates.

Filtering is done using mapped_hgvs_c when it has any non-empty candidate; if
mapped_hgvs_c is absent or entirely blank, mapped_hgvs_g is used instead.
All pipe-delimited list fields aligned to the DNA candidates are compacted in
lockstep so the same candidate positions are removed from every aligned column.

Rows that end up with no retained DNA candidates are still written. If such a
row also has no mapped_hgvs_p value, the script logs a warning count at the end
of execution.
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

logger = logging.getLogger(__name__)

DNA_VARIANT_COLUMNS = [
    "mapped_hgvs_g",
    "mapped_hgvs_c",
    "reverse_translation_warnings",
    "mapped_hgvs_g_chromosome",
    "mapped_hgvs_g_start",
    "mapped_hgvs_g_stop",
    "mapped_hgvs_g_ref",
    "mapped_hgvs_g_alt",
    "mapped_hgvs_c_transcript",
    "mapped_hgvs_c_start",
    "mapped_hgvs_c_stop",
    "mapped_hgvs_c_ref",
    "mapped_hgvs_c_alt",
    "touches_intronic_region",
    "spans_intron",
    "dna_clingen_allele_id",
]

ANNOTATION_PREFIXES = (
    "alphamissense.",
    "clingen_evidence_repository.",
    "clinvar.",
    "gnomad.",
    "mutpred2.",
    "spliceai.",
    "revel.",
    "vep.",
)

_SNV_RE = re.compile(r"^[ACGTN]$", re.IGNORECASE)
_KEEP_CLASS_ALIASES = {
    "wt": {"wt", "wildtype", "wild-type", "no_change", "unchanged"},
    "snv": {"snv", "substitution"},
}


@dataclass(frozen=True)
class FilterStats:
    rows_processed: int
    candidates_removed: int
    rows_without_dna_and_protein: int


def _split_pipe_preserve_positions(value: str) -> list[str]:
    raw = value or ""
    if "|" not in raw:
        return [raw.strip()]
    return [part.strip() for part in raw.split("|")]


def _split_csv_args(values: Sequence[str]) -> list[str]:
    out: list[str] = []
    seen: set[str] = set()
    for raw in values:
        for item in raw.split(","):
            value = item.strip()
            if not value or value in seen:
                continue
            seen.add(value)
            out.append(value)
    return out


def _get_aligned_columns(fieldnames: Sequence[str]) -> list[str]:
    aligned = [col for col in DNA_VARIANT_COLUMNS if col in fieldnames]
    for col in fieldnames:
        if col.startswith(ANNOTATION_PREFIXES) and col not in aligned:
            aligned.append(col)
    return aligned


def _has_non_empty_candidate(values: Sequence[str]) -> bool:
    return any(value.strip() for value in values)


def _normalize_keep_classes(keep_classes_raw: Sequence[str] | None) -> set[str]:
    if not keep_classes_raw:
        return {"snv", "wt"}

    keep_classes: set[str] = set()
    for raw in keep_classes_raw:
        for item in raw.split(","):
            value = item.strip().lower()
            if not value:
                continue
            if value == "all":
                return {"all"}
            alias = next((canonical for canonical, synonyms in _KEEP_CLASS_ALIASES.items() if value in synonyms), value)
            keep_classes.add(alias)
    if not keep_classes:
        return {"snv", "wt"}
    return keep_classes


def _classify_candidate(hgvs: str) -> str:
    value = (hgvs or "").strip()
    if not value:
        return "blank"

    lowered = value.lower()
    if lowered in {"no_change", "no-change"} or value.endswith("="):
        return "wt"

    body = value.split(":", 1)[1] if ":" in value else value
    body = body.strip()

    if "delins" in body:
        return "delins"
    if "dup" in body:
        return "dup"
    if "inv" in body:
        return "inv"
    if "ins" in body:
        return "ins"
    if "del" in body:
        return "del"
    if ">" in body:
        left, right = body.rsplit(">", 1)
        if left and right and _SNV_RE.search(left[-1]) and _SNV_RE.fullmatch(right):
            return "snv"
        return "complex"
    return "other"


def _should_keep_candidate(hgvs: str, keep_classes: set[str]) -> bool:
    if "all" in keep_classes:
        return True
    return _classify_candidate(hgvs) in keep_classes


def filter_dna_variants(
    input_file: Path,
    output_file: Path,
    keep_classes: Sequence[str] | None = None,
    *,
    delimiter: str = "\t",
    skip: int = 0,
    limit: Optional[int] = None,
    csv_field_size_limit: int = csv.field_size_limit(),
) -> FilterStats:
    """Filter DNA candidates in place while preserving the number of rows."""
    if not input_file.exists():
        raise FileNotFoundError(f"Input file not found: {input_file}")
    if skip < 0:
        raise ValueError(f"--skip must be >= 0, got: {skip}")
    if limit is not None and limit < 1:
        raise ValueError(f"--limit must be >= 1 when provided, got: {limit}")
    if not delimiter:
        raise ValueError("--delimiter must not be blank")

    csv.field_size_limit(csv_field_size_limit)
    keep_class_set = _normalize_keep_classes(keep_classes)

    with input_file.open("r", encoding="utf-8", newline="") as in_fh, output_file.open(
        "w", encoding="utf-8", newline=""
    ) as out_fh:
        reader = csv.DictReader(in_fh, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"Input file appears empty: {input_file}")

        fieldnames = list(reader.fieldnames)
        aligned_columns = _get_aligned_columns(fieldnames)
        if not aligned_columns:
            raise ValueError("No DNA variant columns found in input file")

        writer = csv.DictWriter(
            out_fh,
            fieldnames=fieldnames,
            delimiter=delimiter,
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()

        rows_processed = 0
        candidates_removed = 0
        rows_without_dna_and_protein = 0

        selected_rows = reader
        for _ in range(skip):
            next(selected_rows, None)

        for row in selected_rows:
            if limit is not None and rows_processed >= limit:
                break

            c_values = _split_pipe_preserve_positions(row.get("mapped_hgvs_c", "")) if "mapped_hgvs_c" in row else []
            g_values = _split_pipe_preserve_positions(row.get("mapped_hgvs_g", "")) if "mapped_hgvs_g" in row else []
            use_c = _has_non_empty_candidate(c_values)
            primary_values = c_values if use_c else g_values
            candidate_count = max(
                len(primary_values),
                *(len(_split_pipe_preserve_positions(row.get(col, ""))) for col in aligned_columns),
            ) if aligned_columns else len(primary_values)

            if candidate_count == 0:
                candidate_count = 1

            candidate_values_by_column = {
                col: _split_pipe_preserve_positions(row.get(col, "")) for col in aligned_columns
            }

            retained_indices: list[int] = []
            for idx in range(candidate_count):
                source_values = primary_values if primary_values else [""]
                hgvs_value = source_values[idx] if idx < len(source_values) else ""
                if _should_keep_candidate(hgvs_value, keep_class_set):
                    retained_indices.append(idx)

            for col in aligned_columns:
                values = candidate_values_by_column[col]
                retained_values = [values[idx] if idx < len(values) else "" for idx in retained_indices]
                row[col] = "|".join(retained_values) if retained_values else ""

            rows_processed += 1
            candidates_removed += max(0, candidate_count - len(retained_indices))

            has_any_dna = any(
                (row.get(col) or "").strip()
                for col in ("mapped_hgvs_c", "mapped_hgvs_g", "dna_clingen_allele_id")
                if col in row
            )
            if not has_any_dna and not (row.get("mapped_hgvs_p") or "").strip():
                rows_without_dna_and_protein += 1

            writer.writerow(row)

    return FilterStats(
        rows_processed=rows_processed,
        candidates_removed=candidates_removed,
        rows_without_dna_and_protein=rows_without_dna_and_protein,
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("input_file", type=Path, help="Input TSV/CSV file")
    parser.add_argument("output_file", type=Path, help="Output TSV/CSV file")
    parser.add_argument(
        "--keep-class",
        dest="keep_classes_raw",
        action="append",
        default=None,
        metavar="CLASS",
        help=(
            "Retain only candidates in the given class. May be repeated and/or comma-separated. "
            "Supported classes: snv, wt, del, ins, dup, inv, delins, complex, other, all. "
            "Default: snv,wt"
        ),
    )
    parser.add_argument(
        "--delimiter",
        default="\t",
        help="Input/output delimiter (default TAB)",
    )
    parser.add_argument(
        "--skip",
        type=int,
        default=0,
        metavar="N",
        help="Number of data rows to skip before filtering.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        metavar="N",
        help="Maximum number of data rows to process after applying --skip.",
    )
    parser.add_argument(
        "--csv-field-size-limit",
        type=int,
        default=csv.field_size_limit(),
        metavar="BYTES",
        help="Maximum per-field character length for CSV/TSV parsing.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity",
    )
    return parser


def main(argv: Optional[list[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    logging.basicConfig(level=getattr(logging, args.log_level), format="%(asctime)s %(levelname)s %(name)s: %(message)s")

    if args.skip < 0:
        logger.error("--skip must be >= 0, got: %d", args.skip)
        return 1
    if args.limit is not None and args.limit < 1:
        logger.error("--limit must be >= 1 when provided, got: %d", args.limit)
        return 1

    try:
        stats = filter_dna_variants(
            args.input_file,
            args.output_file,
            args.keep_classes_raw,
            delimiter=args.delimiter,
            skip=args.skip,
            limit=args.limit,
            csv_field_size_limit=args.csv_field_size_limit,
        )
    except (FileNotFoundError, ValueError) as exc:
        logger.error(str(exc))
        return 1

    if stats.rows_without_dna_and_protein:
        logger.warning(
            "%d row(s) ended with no retained DNA variants and no mapped_hgvs_p",
            stats.rows_without_dna_and_protein,
        )
    logger.info(
        "Processed %d row(s); removed %d candidate(s); wrote %s",
        stats.rows_processed,
        stats.candidates_removed,
        args.output_file,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
