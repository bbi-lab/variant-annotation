"""Unified CLI for lightweight table utilities.

This module exposes multiple utility commands behind a single entrypoint so
scripts and Docker wrappers can invoke:

    python -m src.utilities <command> [options]

Available commands:
- compare-columns: find rows where paired columns contain differing values
- filter-columns: keep or omit selected columns from CSV/TSV files
- filter-rows: keep rows where selected columns contain values
- merge-rows: merge rows from multiple CSV/TSV files by composite key
"""

from __future__ import annotations

from typing import Optional

import click

from src.compare_columns import compare_columns
from src.filter_columns import filter_columns
from src.filter_rows import filter_rows
from src.merge_rows import merge_rows


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


@click.group()
def main() -> None:
    """Run utility subcommands.

    The first positional argument selects the command.
    """


@main.command("filter-columns")
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
def filter_columns_cmd(
    input_file: str,
    output_file: str,
    keep_cols_raw: tuple[str, ...],
    omit_cols_raw: tuple[str, ...],
) -> None:
    """Filter columns in INPUT_FILE and write OUTPUT_FILE."""
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


@main.command("merge-rows")
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
def merge_rows_cmd(
    output_file: str,
    input_files: tuple[str, ...],
    key_cols_raw: tuple[str, ...],
) -> None:
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


@main.command("compare-columns")
@click.argument("input_file")
@click.option(
    "--col-a",
    "cols_a_raw",
    multiple=True,
    required=True,
    metavar="COLUMN",
    help="Column name for group A. May be repeated to compare multiple pairs.",
)
@click.option(
    "--col-b",
    "cols_b_raw",
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
def compare_columns_cmd(
    input_file: str,
    cols_a_raw: tuple[str, ...],
    cols_b_raw: tuple[str, ...],
    output_file: Optional[str],
    skip: int,
    limit: Optional[int],
) -> None:
    """Find rows where paired column values differ.

    Compares each --col-a column against the corresponding --col-b column
    (by position). Rows with at least one differing pair are written to
    --output (or stdout) with a "differences" column listing the mismatched pairs.
    """
    if len(cols_a_raw) != len(cols_b_raw):
        raise click.ClickException(
            f"--col-a and --col-b must have the same number of entries "
            f"({len(cols_a_raw)} vs {len(cols_b_raw)})."
        )
    try:
        n_diffs = compare_columns(
            input_file=input_file,
            cols_a=list(cols_a_raw),
            cols_b=list(cols_b_raw),
            output_file=output_file,
            skip=skip,
            limit=limit,
        )
        if output_file:
            click.echo(f"Wrote {n_diffs} differing row(s) to {output_file}")
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc


@main.command("filter-rows")
@click.argument("input_file")
@click.argument("output_file")
@click.option(
    "--value-col",
    "value_cols_raw",
    multiple=True,
    required=True,
    metavar="COLUMN",
    help=(
        "Column to test for non-empty value. May be repeated and/or comma-separated, "
        "for example: --value-col a --value-col b or --value-col a,b"
    ),
)
@click.option(
    "--match",
    "match_mode",
    type=click.Choice(["any", "all"], case_sensitive=False),
    default="any",
    show_default=True,
    help="When multiple --value-col columns are supplied, keep rows matching any or all.",
)
def filter_rows_cmd(
    input_file: str,
    output_file: str,
    value_cols_raw: tuple[str, ...],
    match_mode: str,
) -> None:
    """Filter rows in INPUT_FILE and write OUTPUT_FILE."""
    value_columns = _split_csv_args(value_cols_raw)
    match_mode = match_mode.lower()
    try:
        n_rows = filter_rows(
            input_file=input_file,
            output_file=output_file,
            value_columns=value_columns,
            match=match_mode,
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Wrote {n_rows} row(s) to {output_file}")


if __name__ == "__main__":
    main()
