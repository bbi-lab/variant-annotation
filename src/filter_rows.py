"""Filter rows by value presence/absence in selected columns.

Usage::

    python -m src.filter_rows in.tsv out.tsv --value-col col_a
    python -m src.filter_rows in.tsv out.tsv --value-col col_a,col_b --match all
    python -m src.filter_rows in.tsv out.tsv --value-col mapped_hgvs_p --value-state blank
"""

import csv
from pathlib import Path

import click


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


def _has_value(value: str | None) -> bool:
    """Return True if value is non-empty after trimming whitespace."""
    return bool((value or "").strip())


def filter_rows(
    input_file: str,
    output_file: str,
    value_columns: list[str],
    match: str = "any",
    value_state: str = "non-empty",
) -> int:
    """Filter rows by value presence/absence in one or more columns.

    Args:
        input_file: Input CSV/TSV file.
        output_file: Output CSV/TSV file.
        value_columns: Columns to test.
        match: "any" keeps rows where at least one column has a value,
            "all" keeps rows where every specified column has a value.
        value_state: "non-empty" keeps rows where selected columns have values;
            "blank" keeps rows where selected columns are blank.

    Returns:
        Number of rows written.

    Raises:
        ValueError: On invalid input, unknown columns, or bad match mode.
    """
    if not value_columns:
        raise ValueError("At least one --value-col is required.")
    if match not in {"any", "all"}:
        raise ValueError("match must be 'any' or 'all'.")
    if value_state not in {"non-empty", "blank"}:
        raise ValueError("value_state must be 'non-empty' or 'blank'.")

    in_sep = _detect_separator(input_file)
    out_sep = _detect_separator(output_file)

    with open(input_file, newline="") as in_fh:
        reader = csv.DictReader(in_fh, delimiter=in_sep)
        fieldnames = list(reader.fieldnames or [])
        if not fieldnames:
            raise ValueError(f"Input file {input_file!r} is empty or missing a header.")

        missing = [c for c in value_columns if c not in fieldnames]
        if missing:
            raise ValueError(f"Requested value columns not found: {missing}")

        with open(output_file, "w", newline="") as out_fh:
            writer = csv.DictWriter(
                out_fh,
                fieldnames=fieldnames,
                delimiter=out_sep,
                extrasaction="ignore",
            )
            writer.writeheader()

            n_rows = 0
            for row in reader:
                checks = [_has_value(row.get(col)) for col in value_columns]
                if value_state == "blank":
                    checks = [not c for c in checks]
                keep = any(checks) if match == "any" else all(checks)
                if keep:
                    writer.writerow(row)
                    n_rows += 1

    return n_rows


@click.command()
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
def main(
    input_file: str,
    output_file: str,
    value_cols_raw: tuple[str, ...],
    match_mode: str,
    value_state: str,
) -> None:
    """Filter rows in INPUT_FILE and write OUTPUT_FILE."""
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
