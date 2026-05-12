"""Unified CLI for lightweight table utilities.

This module exposes multiple utility commands behind a single entrypoint so
scripts and Docker wrappers can invoke:

    python -m src.utilities <command> [options]

Available commands:
- compare-columns: find rows where paired columns contain differing values
- filter-columns: keep or omit selected columns from CSV/TSV files
- filter-rows: keep rows where selected columns contain values
- replace-rows: replace rows from multiple CSV/TSV files by composite key
- merge-columns: left-join two files and add selected columns from the second
- reorder-columns: reorder columns using a requested column order
"""

from __future__ import annotations

import csv
from typing import Optional

import click

from src.compare_columns import compare_columns
from src.filter_columns import filter_columns, rename_columns, _parse_keep_col_args
from src.filter_rows import filter_rows
from src.merge_columns import merge_columns, _parse_add_col_args, _parse_key_col_args
from src.reorder_columns import reorder_columns
from src.replace_rows import replace_rows


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
@click.option(
    "--csv-field-size-limit",
    default=csv.field_size_limit(),
    show_default=True,
    type=click.IntRange(min=1),
    help="Maximum per-field character length for CSV/TSV parsing.",
)
def filter_columns_cmd(
    input_file: str,
    output_file: str,
    keep_cols_raw: tuple[str, ...],
    omit_cols_raw: tuple[str, ...],
    csv_field_size_limit: int,
) -> None:
    """Filter columns in INPUT_FILE and write OUTPUT_FILE."""
    csv.field_size_limit(csv_field_size_limit)
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


@main.command("replace-rows")
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
@click.option(
    "--csv-field-size-limit",
    default=csv.field_size_limit(),
    show_default=True,
    type=click.IntRange(min=1),
    help="Maximum per-field character length for CSV/TSV parsing.",
)
def replace_rows_cmd(
    output_file: str,
    input_files: tuple[str, ...],
    key_cols_raw: tuple[str, ...],
    csv_field_size_limit: int,
) -> None:
    """Replace rows from INPUT_FILES into OUTPUT_FILE using composite keys.

    Later files override matching keys from earlier files. New keys from later
    files are appended.
    """
    csv.field_size_limit(csv_field_size_limit)
    key_columns = _split_csv_args(key_cols_raw)
    try:
        n_rows = replace_rows(
            input_files=list(input_files),
            output_file=output_file,
            key_columns=key_columns,
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Wrote {n_rows} row(s) to {output_file}")


@main.command("merge-columns")
@click.argument("base_file")
@click.argument("extra_file")
@click.argument("output_file")
@click.option(
    "--key-col",
    "key_cols_raw",
    multiple=True,
    required=True,
    metavar="BASE_COL[:EXTRA_COL]",
    help=(
        "Key column(s) for join identity. May be repeated and/or comma-separated. "
        "Use BASE_COL:EXTRA_COL when the column has different names in the two files "
        "(e.g. --key-col 'dataset_name:Dataset Name'). "
        "Plain names mean the same column is used on both sides."
    ),
)
@click.option(
    "--add-col",
    "add_cols_raw",
    multiple=True,
    metavar="COLUMN[:OUTPUT_NAME]",
    help=(
        "Column(s) to add from extra file. May be repeated and/or comma-separated. "
        "Use SRC:DEST to rename a column in the output "
        "(e.g. --add-col 'Detects Splicing Variants?:splice_measure'). "
        "If omitted, all new non-key columns from extra are added."
    ),
)
@click.option(
    "--add-all-cols-from-extra",
    is_flag=True,
    help="Add all non-key columns from extra that are not already in base.",
)
@click.option(
    "--extra-worksheet",
    default=None,
    metavar="NAME_OR_INDEX",
    help=(
        "Worksheet to read when EXTRA_FILE is an Excel workbook (.xlsx/.xls). "
        "Accepts a sheet name or a 1-based integer index. "
        "Defaults to the first (active) sheet."
    ),
)
@click.option(
    "--csv-field-size-limit",
    default=csv.field_size_limit(),
    show_default=True,
    type=click.IntRange(min=1),
    help="Maximum per-field character length for CSV/TSV parsing.",
)
def merge_columns_cmd(
    base_file: str,
    extra_file: str,
    output_file: str,
    key_cols_raw: tuple[str, ...],
    add_cols_raw: tuple[str, ...],
    add_all_cols_from_extra: bool,
    extra_worksheet: str | None,
    csv_field_size_limit: int,
) -> None:
    """Left-join BASE_FILE with EXTRA_FILE and write OUTPUT_FILE."""
    csv.field_size_limit(csv_field_size_limit)
    key_columns, extra_key_columns = _parse_key_col_args(key_cols_raw)
    add_columns, column_renames = _parse_add_col_args(add_cols_raw)

    try:
        n_rows = merge_columns(
            base_file=base_file,
            extra_file=extra_file,
            output_file=output_file,
            key_columns=key_columns,
            add_columns=add_columns,
            add_all_cols_from_extra=add_all_cols_from_extra,
            column_renames=column_renames,
            extra_worksheet=extra_worksheet,
            extra_key_columns=extra_key_columns,
        )
    except (ValueError, ImportError) as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Wrote {n_rows} row(s) to {output_file}")


@main.command("rename-columns")
@click.argument("input_file")
@click.argument("output_file")
@click.option(
    "--keep-col",
    "keep_cols_raw",
    multiple=True,
    metavar="SRC[:DEST]",
    help=(
        "Column to keep, optionally renaming it. May be repeated and/or comma-separated. "
        "Use SRC:DEST to rename (e.g. --keep-col 'old_name:new_name'). "
        "Mutually exclusive with --omit-col."
    ),
)
@click.option(
    "--omit-col",
    "omit_cols_raw",
    multiple=True,
    metavar="COLUMN",
    help=(
        "Column to omit, keeping all others unchanged. May be repeated and/or "
        "comma-separated. Mutually exclusive with --keep-col."
    ),
)
@click.option(
    "--reorder",
    is_flag=True,
    help=(
        "Output columns in the order given by --keep-col rather than input file order. "
        "Only applies in keep mode."
    ),
)
@click.option(
    "--csv-field-size-limit",
    default=csv.field_size_limit(),
    show_default=True,
    type=click.IntRange(min=1),
    help="Maximum per-field character length for CSV/TSV parsing.",
)
def rename_columns_cmd(
    input_file: str,
    output_file: str,
    keep_cols_raw: tuple[str, ...],
    omit_cols_raw: tuple[str, ...],
    reorder: bool,
    csv_field_size_limit: int,
) -> None:
    """Keep/omit and optionally rename columns in INPUT_FILE, writing OUTPUT_FILE.

    In keep mode (--keep-col), only the listed columns are written.  Use
    SRC:DEST syntax to rename a column at the same time.  Pass --reorder to
    output columns in the order given by --keep-col.

    In omit mode (--omit-col), all columns except the listed ones are written,
    names unchanged.
    """
    csv.field_size_limit(csv_field_size_limit)
    keep_specs = _parse_keep_col_args(keep_cols_raw)
    omit_cols = _split_csv_args(omit_cols_raw)
    try:
        n_rows = rename_columns(
            input_file=input_file,
            output_file=output_file,
            keep_specs=keep_specs,
            omit_cols=omit_cols,
            reorder=reorder,
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Wrote {n_rows} row(s) to {output_file}")


@main.command("reorder-columns")
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
@click.option(
    "--csv-field-size-limit",
    default=csv.field_size_limit(),
    show_default=True,
    type=click.IntRange(min=1),
    help="Maximum per-field character length for CSV/TSV parsing.",
)
def reorder_columns_cmd(
    input_file: str,
    output_file: str,
    column_order_raw: tuple[str, ...],
    csv_field_size_limit: int,
) -> None:
    """Reorder columns in INPUT_FILE and write OUTPUT_FILE."""
    csv.field_size_limit(csv_field_size_limit)
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
@click.option(
    "--csv-field-size-limit",
    default=csv.field_size_limit(),
    show_default=True,
    type=click.IntRange(min=1),
    help="Maximum per-field character length for CSV/TSV parsing.",
)
def compare_columns_cmd(
    input_file: str,
    cols_a_raw: tuple[str, ...],
    cols_b_raw: tuple[str, ...],
    output_file: Optional[str],
    skip: int,
    limit: Optional[int],
    csv_field_size_limit: int,
) -> None:
    """Find rows where paired column values differ.

    Compares each --col-a column against the corresponding --col-b column
    (by position). Rows with at least one differing pair are written to
    --output (or stdout) with a "differences" column listing the mismatched pairs.
    """
    csv.field_size_limit(csv_field_size_limit)
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
@click.option(
    "--value-state",
    type=click.Choice(["non-empty", "blank"], case_sensitive=False),
    default="non-empty",
    show_default=True,
    help="Whether selected columns should be non-empty or blank.",
)
@click.option(
    "--csv-field-size-limit",
    default=csv.field_size_limit(),
    show_default=True,
    type=click.IntRange(min=1),
    help="Maximum per-field character length for CSV/TSV parsing.",
)
def filter_rows_cmd(
    input_file: str,
    output_file: str,
    value_cols_raw: tuple[str, ...],
    match_mode: str,
    value_state: str,
    csv_field_size_limit: int,
) -> None:
    """Filter rows in INPUT_FILE and write OUTPUT_FILE."""
    csv.field_size_limit(csv_field_size_limit)
    value_columns = _split_csv_args(value_cols_raw)
    match_mode = match_mode.lower()
    value_state = value_state.lower()
    try:
        n_rows = filter_rows(
            input_file=input_file,
            output_file=output_file,
            value_columns=value_columns,
            match=match_mode,
            value_state=value_state,
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Wrote {n_rows} row(s) to {output_file}")


if __name__ == "__main__":
    main()
