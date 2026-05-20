"""
Compare the mapped_hgvs_g, mapped_hgvs_c, and mapped_hgvs_p columns between
two TSV files, joining on variant_urn.

For each pipe-separated component that differs, outputs one row with:
  variant_urn, column, index, old_value, new_value

If a variant_urn exists in only one file, the missing side's component values
are blank.

Output defaults to <new_tsv_stem>.hgvs_diffs.tsv in the current directory.
Override with -o / --output.

Both files must contain a variant_urn column.
Runs locally using the project virtual environment (not Docker).
"""

import argparse
import csv
import sys
from pathlib import Path

csv.field_size_limit(10 * 1024 * 1024)  # 10 MB

HGVS_COLS = ("mapped_hgvs_g", "mapped_hgvs_c", "mapped_hgvs_p")
KEY_COL = "variant_urn"
DIFF_HEADERS = ["variant_urn", "column", "index", "old_value", "new_value"]


def load_hgvs(path: str) -> dict[str, dict[str, str]]:
    """Return {variant_urn: {col: value}} for the three HGVS columns."""
    rows: dict[str, dict[str, str]] = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Empty file: {path}")
        for row in reader:
            key = row[KEY_COL]
            rows[key] = {col: row.get(col, "") for col in HGVS_COLS}
    return rows


def component_diffs(
    urn: str,
    col: str,
    old_val: str,
    new_val: str,
) -> list[dict[str, str]]:
    """Compare pipe-separated components; return a row for each that differs."""
    old_parts = old_val.split("|") if old_val else []
    new_parts = new_val.split("|") if new_val else []
    length = max(len(old_parts), len(new_parts))
    diffs = []
    for i in range(length):
        old_comp = old_parts[i] if i < len(old_parts) else ""
        new_comp = new_parts[i] if i < len(new_parts) else ""
        if old_comp != new_comp:
            diffs.append(
                {
                    "variant_urn": urn,
                    "column": col,
                    "index": str(i),
                    "old_value": old_comp,
                    "new_value": new_comp,
                }
            )
    return diffs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("old_tsv", help="Path to the old (baseline) TSV file")
    parser.add_argument("new_tsv", help="Path to the new TSV file")
    parser.add_argument(
        "-o",
        "--output",
        help="Output file path (default: <new_tsv_stem>.hgvs_diffs.tsv)",
    )
    args = parser.parse_args()

    output_path = args.output or f"{Path(args.new_tsv).stem}.hgvs_diffs.tsv"

    print(f"Loading {args.old_tsv} ...", file=sys.stderr)
    old_rows = load_hgvs(args.old_tsv)
    print(f"Loading {args.new_tsv} ...", file=sys.stderr)
    new_rows = load_hgvs(args.new_tsv)

    all_urns = old_rows.keys() | new_rows.keys()
    empty = {col: "" for col in HGVS_COLS}

    diff_rows: list[dict[str, str]] = []
    for urn in sorted(all_urns):
        old = old_rows.get(urn, empty)
        new = new_rows.get(urn, empty)
        for col in HGVS_COLS:
            if old[col] != new[col]:
                diff_rows.extend(component_diffs(urn, col, old[col], new[col]))

    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=DIFF_HEADERS, delimiter="\t")
        writer.writeheader()
        writer.writerows(diff_rows)

    print(f"Wrote {len(diff_rows)} HGVS component diffs → {output_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
