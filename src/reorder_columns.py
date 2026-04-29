"""Reorder CSV/TSV columns.

Usage::

    python -m src.reorder_columns in.tsv out.tsv --column-order id,gene,value
"""

from __future__ import annotations

import csv
from pathlib import Path

import click


def _detect_separator(file_path: str) -> str:
    return "\t" if Path(file_path).suffix.lower() in (".tsv", ".txt") else ","


def _split_csv_args(values: tuple[str, ...]) -> list[str]:
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


def reorder_columns(
    input_file: str,
    output_file: str,
    column_order: list[str],
) -> int:
    """Write *input_file* to *output_file* with columns in *column_order*.

    Columns not listed in *column_order* are appended at the end in their
    original order.
    """
    if not column_order:
        raise ValueError("At least one column in --column-order is required.")

    sep = _detect_separator(input_file)
    with open(input_file, newline="", encoding="utf-8") as in_fh:
        reader = csv.DictReader(in_fh, delimiter=sep)
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            raise ValueError(f"Input file {input_file!r} is empty or missing a header.")

        missing = [c for c in column_order if c not in fieldnames]
        if missing:
            raise ValueError(f"Input file is missing columns from --column-order: {missing}")

        out_fieldnames = list(column_order) + [c for c in fieldnames if c not in column_order]
        rows = list(reader)

    out_sep = _detect_separator(output_file)
    with open(output_file, "w", newline="", encoding="utf-8") as out_fh:
        writer = csv.DictWriter(
            out_fh,
            fieldnames=out_fieldnames,
            delimiter=out_sep,
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

    return len(rows)


@click.command()
@click.argument("input_file")
@click.argument("output_file")
@click.option(
    "--column-order",
    "column_order_raw",
    multiple=True,
    required=True,
    metavar="COLUMN",
    help=(
        "Desired column order. May be repeated and/or comma-separated, "
        "for example: --column-order a,b,c"
    ),
)
def main(input_file: str, output_file: str, column_order_raw: tuple[str, ...]) -> None:
    column_order = _split_csv_args(column_order_raw)
    try:
        n_rows = reorder_columns(
            input_file=input_file,
            output_file=output_file,
            column_order=column_order,
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Wrote {n_rows} row(s) to {output_file}")


if __name__ == "__main__":
    main()
