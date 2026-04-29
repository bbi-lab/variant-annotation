"""Merge two CSV/TSV files by key columns, selecting columns from each side.

Typical use is to keep all columns from a base file and add one or more columns
from a second file by joining on key column(s).

Usage::

    python -m src.merge_columns base.tsv extra.tsv out.tsv --key-col id --add-col score
    python -m src.merge_columns base.tsv extra.tsv out.tsv --key-col gene,variant --add-all-cols-from-extra
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


def _build_key(row: dict[str, str], key_columns: list[str]) -> tuple[str, ...]:
    return tuple((row.get(col) or "").strip() for col in key_columns)


def merge_columns(
    base_file: str,
    extra_file: str,
    output_file: str,
    key_columns: list[str],
    add_columns: list[str],
    add_all_cols_from_extra: bool = False,
) -> int:
    """Left-join *base_file* with *extra_file* and write selected merged columns.

    The output always contains all base-file rows. The join uses *key_columns*.
    """
    if not key_columns:
        raise ValueError("At least one key column is required.")

    base_sep = _detect_separator(base_file)
    extra_sep = _detect_separator(extra_file)

    with open(base_file, newline="", encoding="utf-8") as base_fh:
        base_reader = csv.DictReader(base_fh, delimiter=base_sep)
        base_fieldnames = list(base_reader.fieldnames or [])
        if not base_fieldnames:
            raise ValueError(f"Input file {base_file!r} is empty or missing a header.")
        missing_base_keys = [c for c in key_columns if c not in base_fieldnames]
        if missing_base_keys:
            raise ValueError(f"Base file is missing key columns: {missing_base_keys}")
        base_rows = list(base_reader)

    with open(extra_file, newline="", encoding="utf-8") as extra_fh:
        extra_reader = csv.DictReader(extra_fh, delimiter=extra_sep)
        extra_fieldnames = list(extra_reader.fieldnames or [])
        if not extra_fieldnames:
            raise ValueError(f"Input file {extra_file!r} is empty or missing a header.")
        missing_extra_keys = [c for c in key_columns if c not in extra_fieldnames]
        if missing_extra_keys:
            raise ValueError(f"Extra file is missing key columns: {missing_extra_keys}")

        if add_all_cols_from_extra:
            selected_add_columns = [c for c in extra_fieldnames if c not in key_columns and c not in base_fieldnames]
        else:
            if not add_columns:
                selected_add_columns = [c for c in extra_fieldnames if c not in key_columns and c not in base_fieldnames]
            else:
                missing_add_cols = [c for c in add_columns if c not in extra_fieldnames]
                if missing_add_cols:
                    raise ValueError(f"Extra file is missing add columns: {missing_add_cols}")
                selected_add_columns = add_columns

        extra_by_key: dict[tuple[str, ...], dict[str, str]] = {}
        for row in extra_reader:
            extra_by_key[_build_key(row, key_columns)] = row

    output_fieldnames = list(base_fieldnames)
    for col in selected_add_columns:
        if col not in output_fieldnames:
            output_fieldnames.append(col)

    out_sep = _detect_separator(output_file)
    with open(output_file, "w", newline="", encoding="utf-8") as out_fh:
        writer = csv.DictWriter(
            out_fh,
            fieldnames=output_fieldnames,
            delimiter=out_sep,
            extrasaction="ignore",
        )
        writer.writeheader()

        for base_row in base_rows:
            merged = dict(base_row)
            extra_row = extra_by_key.get(_build_key(base_row, key_columns))
            if extra_row is not None:
                for col in selected_add_columns:
                    merged[col] = extra_row.get(col, "")
            else:
                for col in selected_add_columns:
                    merged.setdefault(col, "")
            writer.writerow(merged)

    return len(base_rows)


@click.command()
@click.argument("base_file")
@click.argument("extra_file")
@click.argument("output_file")
@click.option(
    "--key-col",
    "key_cols_raw",
    multiple=True,
    required=True,
    metavar="COLUMN",
    help=(
        "Key column(s) for join identity. May be repeated and/or comma-separated, "
        "for example: --key-col gene --key-col variant or --key-col gene,variant"
    ),
)
@click.option(
    "--add-col",
    "add_cols_raw",
    multiple=True,
    metavar="COLUMN",
    help=(
        "Column(s) to add from extra file. May be repeated and/or comma-separated. "
        "If omitted, all new non-key columns from extra are added."
    ),
)
@click.option(
    "--add-all-cols-from-extra",
    is_flag=True,
    help="Add all non-key columns from extra that are not already in base.",
)
def main(
    base_file: str,
    extra_file: str,
    output_file: str,
    key_cols_raw: tuple[str, ...],
    add_cols_raw: tuple[str, ...],
    add_all_cols_from_extra: bool,
) -> None:
    key_columns = _split_csv_args(key_cols_raw)
    add_columns = _split_csv_args(add_cols_raw)

    try:
        n_rows = merge_columns(
            base_file=base_file,
            extra_file=extra_file,
            output_file=output_file,
            key_columns=key_columns,
            add_columns=add_columns,
            add_all_cols_from_extra=add_all_cols_from_extra,
        )
    except ValueError as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(f"Wrote {n_rows} row(s) to {output_file}")


if __name__ == "__main__":
    main()
