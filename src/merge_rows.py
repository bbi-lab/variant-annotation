"""Merge rows from multiple CSV/TSV files using a composite key.

Rows are processed in file order. For duplicate keys, later rows replace earlier
rows. New keys from later files are appended to the output.

Usage::

    python -m src.merge_rows out.tsv in1.tsv in2.tsv in3.tsv --key-col id
    python -m src.merge_rows out.csv base.csv patch.csv --key-col gene,variant
"""

import csv
from pathlib import Path

import click


# ---------------------------------------------------------------------------
# File / format utilities
# ---------------------------------------------------------------------------


def _detect_separator(file_path: str) -> str:
    """Return field delimiter based on file extension."""
    return "\t" if Path(file_path).suffix.lower() in (".tsv", ".txt") else ","


def _split_csv_args(values: tuple[str, ...]) -> list[str]:
    """Flatten repeated/comma-delimited CLI values into a unique list preserving order."""
    out: list[str] = []
    seen: set[str] = set()
    for raw in values:
        for item in raw.split(","):
            col = item.strip()
            if not col:
                continue
            if col not in seen:
                seen.add(col)
                out.append(col)
    return out


def _build_key(row: dict[str, str], key_columns: list[str]) -> tuple[str, ...]:
    """Build a stable tuple key for a row."""
    return tuple((row.get(col) or "").strip() for col in key_columns)


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------


def merge_rows(
    input_files: list[str],
    output_file: str,
    key_columns: list[str],
) -> int:
    """Merge rows from *input_files* and write to *output_file*.

    Args:
        input_files: Ordered list of input CSV/TSV files.
        output_file: Output CSV/TSV file path.
        key_columns: Columns that jointly define a unique row key.

    Returns:
        Number of merged rows written.

    Raises:
        ValueError: If no input files are provided, key columns are empty,
            key columns are missing, or input columns do not match.
    """
    if not input_files:
        raise ValueError("At least one input file is required.")
    if not key_columns:
        raise ValueError("At least one key column is required.")

    canonical_fieldnames: list[str] | None = None
    merged_by_key: dict[tuple[str, ...], dict[str, str]] = {}

    for file_idx, file_path in enumerate(input_files):
        sep = _detect_separator(file_path)
        with open(file_path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter=sep)
            fieldnames = list(reader.fieldnames or [])
            if not fieldnames:
                raise ValueError(f"Input file {file_path!r} is empty or missing a header.")

            missing_keys = [c for c in key_columns if c not in fieldnames]
            if missing_keys:
                raise ValueError(
                    f"Input file {file_path!r} is missing key columns: {missing_keys}"
                )

            if file_idx == 0:
                canonical_fieldnames = fieldnames
            else:
                assert canonical_fieldnames is not None
                if set(fieldnames) != set(canonical_fieldnames):
                    raise ValueError(
                        f"Input file {file_path!r} columns do not match first file. "
                        f"Expected {canonical_fieldnames}, got {fieldnames}"
                    )

            for row in reader:
                key = _build_key(row, key_columns)
                merged_by_key[key] = row

    assert canonical_fieldnames is not None
    out_sep = _detect_separator(output_file)
    with open(output_file, "w", newline="") as out_fh:
        writer = csv.DictWriter(
            out_fh,
            fieldnames=canonical_fieldnames,
            delimiter=out_sep,
            extrasaction="ignore",
        )
        writer.writeheader()
        for row in merged_by_key.values():
            writer.writerow(row)

    return len(merged_by_key)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


@click.command()
@click.argument("output_file")
@click.argument("input_files", nargs=-1, required=True)
@click.option(
    "--key-col",
    "key_cols_raw",
    multiple=True,
    required=True,
    metavar="COLUMN",
    help=(
        "Key column for row identity. May be repeated and/or comma-separated, "
        "for example: --key-col gene --key-col variant or --key-col gene,variant"
    ),
)
def main(output_file: str, input_files: tuple[str, ...], key_cols_raw: tuple[str, ...]) -> None:
    """Merge rows from INPUT_FILES into OUTPUT_FILE using composite keys.

    Later files override matching keys from earlier files. New keys from later
    files are appended.
    """
    key_columns = _split_csv_args(key_cols_raw)
    try:
        n_rows = merge_rows(
            input_files=list(input_files),
            output_file=output_file,
            key_columns=key_columns,
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Wrote {n_rows} merged row(s) to {output_file}")


if __name__ == "__main__":
    main()
