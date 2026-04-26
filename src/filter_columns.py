"""Filter CSV/TSV columns by keep-list or omit-list.

Usage::

    python -m src.filter_columns input.tsv output.tsv --keep-col col1 --keep-col col2
    python -m src.filter_columns input.tsv output.tsv --omit-col big_blob --omit-col notes

Exactly one mode must be selected:
- keep mode: include only selected columns (in input order)
- omit mode: include all columns except selected columns
"""

import csv
from pathlib import Path
from typing import Optional

import click


# ---------------------------------------------------------------------------
# File / format utilities
# ---------------------------------------------------------------------------


def _detect_separator(file_path: str) -> str:
    """Return the field delimiter appropriate for *file_path* based on suffix."""
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


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------


def filter_columns(
    input_file: str,
    output_file: str,
    keep_cols: Optional[list[str]] = None,
    omit_cols: Optional[list[str]] = None,
) -> int:
    """Filter columns from *input_file* into *output_file*.

    Exactly one of *keep_cols* or *omit_cols* must be provided.

    Returns:
        Number of data rows written.

    Raises:
        ValueError: If modes are invalid, or requested columns are missing.
    """
    keep_cols = keep_cols or []
    omit_cols = omit_cols or []

    if bool(keep_cols) == bool(omit_cols):
        raise ValueError("Specify exactly one mode: --keep-col (one or more) OR --omit-col (one or more).")

    in_sep = _detect_separator(input_file)
    out_sep = _detect_separator(output_file)

    with open(input_file, newline="") as in_fh:
        reader = csv.DictReader(in_fh, delimiter=in_sep)
        input_fieldnames = list(reader.fieldnames or [])
        if not input_fieldnames:
            raise ValueError(f"Input file {input_file!r} is empty or missing a header.")

        if keep_cols:
            missing = [c for c in keep_cols if c not in input_fieldnames]
            if missing:
                raise ValueError(f"Requested keep columns not found: {missing}")
            output_fieldnames = [c for c in input_fieldnames if c in set(keep_cols)]
        else:
            missing = [c for c in omit_cols if c not in input_fieldnames]
            if missing:
                raise ValueError(f"Requested omit columns not found: {missing}")
            omit_set = set(omit_cols)
            output_fieldnames = [c for c in input_fieldnames if c not in omit_set]
            if not output_fieldnames:
                raise ValueError("Omit mode would remove all columns; at least one output column is required.")

        with open(output_file, "w", newline="") as out_fh:
            writer = csv.DictWriter(
                out_fh,
                fieldnames=output_fieldnames,
                delimiter=out_sep,
                extrasaction="ignore",
            )
            writer.writeheader()

            n_rows = 0
            for row in reader:
                writer.writerow(row)
                n_rows += 1

    return n_rows


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


@click.command()
@click.argument("input_file")
@click.argument("output_file")
@click.option(
    "--keep-col",
    "keep_cols_raw",
    multiple=True,
    metavar="COLUMN",
    help=(
        "Column to keep. May be repeated. Also accepts comma-separated values, "
        "for example: --keep-col a,b,c"
    ),
)
@click.option(
    "--omit-col",
    "omit_cols_raw",
    multiple=True,
    metavar="COLUMN",
    help=(
        "Column to omit. May be repeated. Also accepts comma-separated values, "
        "for example: --omit-col x,y"
    ),
)
def main(
    input_file: str,
    output_file: str,
    keep_cols_raw: tuple[str, ...],
    omit_cols_raw: tuple[str, ...],
) -> None:
    """Filter columns in INPUT_FILE and write OUTPUT_FILE.

    
    Keep mode example:
        python -m src.filter_columns in.tsv out.tsv --keep-col a --keep-col b

    
    Omit mode example:
        python -m src.filter_columns in.tsv out.tsv --omit-col large_blob,notes
    """
    keep_cols = _split_csv_args(keep_cols_raw)
    omit_cols = _split_csv_args(omit_cols_raw)
    try:
        n_rows = filter_columns(
            input_file=input_file,
            output_file=output_file,
            keep_cols=keep_cols,
            omit_cols=omit_cols,
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Wrote {n_rows} row(s) to {output_file}")


if __name__ == "__main__":
    main()
