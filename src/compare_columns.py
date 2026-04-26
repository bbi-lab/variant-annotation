"""Compare paired columns in a CSV/TSV file to find rows where values differ.

This script is intended to validate newly computed mappings against preexisting
reference columns. Given two lists of column names (A and B) of the same length,
it reports every row where at least one A/B pair has differing values.

Usage::

    python -m src.compare_columns INPUT_FILE \\
        --col-a mapped_hgvs_g --col-b ref_hgvs_g \\
        --col-a mapped_hgvs_c --col-b ref_hgvs_c

    # Write differing rows to a file instead of stdout:
    python -m src.compare_columns INPUT_FILE \\
        --col-a mapped_hgvs_g --col-b ref_hgvs_g \\
        --output differences.tsv

Library usage::

    from src.compare_columns import compare_columns

    n_diffs = compare_columns(
        input_file="annotated.tsv",
        cols_a=["mapped_hgvs_g", "mapped_hgvs_c"],
        cols_b=["ref_hgvs_g", "ref_hgvs_c"],
        output_file="differences.tsv",
    )
"""

import csv
import sys
from pathlib import Path
from typing import Optional

import click


# ---------------------------------------------------------------------------
# File / format utilities
# ---------------------------------------------------------------------------


def _detect_separator(file_path: str) -> str:
    """Return the field delimiter appropriate for *file_path* based on its suffix."""
    return "\t" if Path(file_path).suffix.lower() in (".tsv", ".txt") else ","


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------


def compare_columns(
    input_file: str,
    cols_a: list[str],
    cols_b: list[str],
    output_file: Optional[str] = None,
    skip: int = 0,
    limit: Optional[int] = None,
) -> int:
    """Find rows where any A/B column pair contains differing values.

    Args:
        input_file: Path to the input CSV or TSV file.
        cols_a: List of column names to treat as group A.
        cols_b: List of column names to treat as group B (paired with *cols_a*).
        output_file: If given, write differing rows here as TSV (with a
            ``differences`` column). If None, write to stdout.
        skip: Number of data rows to skip from the start of the file.
        limit: Maximum number of data rows to examine.

    Returns:
        The number of rows that had at least one differing pair.

    Raises:
        ValueError: If *cols_a* and *cols_b* have different lengths, or if any
            column is missing from the input file.
    """
    if len(cols_a) != len(cols_b):
        raise ValueError(
            f"--col-a and --col-b must have the same number of entries "
            f"({len(cols_a)} vs {len(cols_b)})."
        )

    in_sep = _detect_separator(input_file)

    with open(input_file, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=in_sep)
        all_rows = list(reader)
        fieldnames: list[str] = list(reader.fieldnames or [])

    # Validate that all requested columns exist.
    missing = [c for c in cols_a + cols_b if c not in fieldnames]
    if missing:
        raise ValueError(
            f"The following columns were not found in {input_file!r}: {missing}"
        )

    rows = all_rows[skip : (skip + limit) if limit is not None else None]

    # Build output field list: all original columns + a "differences" column.
    out_fieldnames = fieldnames + ["differences"]

    pairs = list(zip(cols_a, cols_b))

    # Choose output destination.
    if output_file:
        out_sep = _detect_separator(output_file)
        out_fh = open(output_file, "w", newline="")
    else:
        out_sep = "\t"
        out_fh = sys.stdout

    n_diffs = 0
    try:
        writer = csv.DictWriter(
            out_fh,
            fieldnames=out_fieldnames,
            delimiter=out_sep,
            extrasaction="ignore",
        )
        writer.writeheader()

        for row_idx, row in enumerate(rows):
            diff_pairs: list[str] = []
            for col_a, col_b in pairs:
                val_a = (row.get(col_a) or "").strip()
                val_b = (row.get(col_b) or "").strip()
                if val_a != val_b:
                    diff_pairs.append(f"{col_a}≠{col_b}")

            if diff_pairs:
                n_diffs += 1
                out_row = dict(row)
                out_row["differences"] = "; ".join(diff_pairs)
                writer.writerow(out_row)
    finally:
        if output_file and out_fh is not sys.stdout:
            out_fh.close()

    # Summary to stderr.
    total = len(rows)
    print(
        f"Found {n_diffs}/{total} rows with differences across {len(pairs)} column pair(s).",
        file=sys.stderr,
    )
    return n_diffs


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


@click.command()
@click.argument("input_file")
@click.option(
    "--col-a",
    "cols_a",
    multiple=True,
    required=True,
    metavar="COLUMN",
    help="Column name for group A. May be repeated to compare multiple pairs.",
)
@click.option(
    "--col-b",
    "cols_b",
    multiple=True,
    required=True,
    metavar="COLUMN",
    help="Column name for group B. Must be repeated the same number of times as --col-a.",
)
@click.option(
    "--output",
    "output_file",
    default=None,
    metavar="FILE",
    help="Write differing rows to FILE instead of stdout.",
)
@click.option(
    "--skip",
    default=0,
    show_default=True,
    type=int,
    help="Number of data rows to skip from the start of the file.",
)
@click.option(
    "--limit",
    default=None,
    type=int,
    help="Maximum number of data rows to examine. Examines all rows when omitted.",
)
def main(
    input_file: str,
    cols_a: tuple[str, ...],
    cols_b: tuple[str, ...],
    output_file: Optional[str],
    skip: int,
    limit: Optional[int],
) -> None:
    """Find rows in INPUT_FILE where paired column values differ.

    Compares each --col-a column against the corresponding --col-b column
    (by position). Rows with at least one differing pair are written to
    --output (or stdout) as a TSV with a "differences" column listing the
    mismatched pairs.

    \b
    Example – compare new and old HGVS mappings in a single file:
        compare_columns annotated.tsv \\
            --col-a mapped_hgvs_g --col-b old_hgvs_g \\
            --col-a mapped_hgvs_c --col-b old_hgvs_c \\
            --output diffs.tsv
    """
    try:
        compare_columns(
            input_file=input_file,
            cols_a=list(cols_a),
            cols_b=list(cols_b),
            output_file=output_file,
            skip=skip,
            limit=limit,
        )
    except ValueError as exc:
        raise click.UsageError(str(exc)) from exc


if __name__ == "__main__":
    main()
