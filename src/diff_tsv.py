"""
Diff two TSV files, comparing rows by variant_urn.

Both files must contain a variant_urn column, which is used as the join key.
Only rows present in the new file are considered; rows that exist only in the
old file are not reported.

Outputs:
  (a) <output_prefix>.changed_rows.tsv
        All rows from the new file that are either new (no matching variant_urn
        in the old file) or differ in any column value from the old file.
        Contains all columns from the new file, in their original order.

  (b) <output_prefix>.hgvs_diffs.tsv
        Per-component diffs of mapped_hgvs_g, mapped_hgvs_c, and mapped_hgvs_p
        for rows present in both files. One row per changed pipe-separated
        component: variant_urn, column, index, old_value, new_value.
        Only produced for rows whose HGVS columns actually changed.

The output prefix defaults to the stem of the new TSV filename.
Override with -o / --output-prefix.

Runs locally using the project virtual environment (not Docker).
"""

import argparse
import csv
import sys
from pathlib import Path

# Some HGVS fields are very large; raise the field-size limit accordingly.
csv.field_size_limit(10 * 1024 * 1024)  # 10 MB

HGVS_COLS = ("mapped_hgvs_g", "mapped_hgvs_c", "mapped_hgvs_p")
KEY_COL = "variant_urn"


def load_tsv(path: str) -> tuple[list[str], dict[str, dict[str, str]]]:
    """Return (headers, {variant_urn: row_dict})."""
    rows: dict[str, dict[str, str]] = {}
    headers: list[str] = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Empty file: {path}")
        headers = list(reader.fieldnames)
        for row in reader:
            key = row[KEY_COL]
            rows[key] = dict(row)
    return headers, rows


def component_diffs(
    urn: str,
    col: str,
    old_val: str,
    new_val: str,
) -> list[dict[str, str]]:
    """Compare pipe-separated components; return rows for any that differ."""
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
        "--output-prefix",
        help="Prefix for output files (default: derived from new_tsv)",
    )
    args = parser.parse_args()

    prefix = args.output_prefix or Path(args.new_tsv).stem
    changed_path = f"{prefix}.changed_rows.tsv"
    diffs_path = f"{prefix}.hgvs_diffs.tsv"

    print(f"Loading {args.old_tsv} ...", file=sys.stderr)
    old_headers, old_rows = load_tsv(args.old_tsv)
    print(f"Loading {args.new_tsv} ...", file=sys.stderr)
    new_headers, new_rows = load_tsv(args.new_tsv)

    changed_rows: list[dict[str, str]] = []
    hgvs_diff_rows: list[dict[str, str]] = []

    for urn, new_row in new_rows.items():
        old_row = old_rows.get(urn)
        if old_row is None or new_row != old_row:
            changed_rows.append(new_row)
            if old_row is not None:
                for col in HGVS_COLS:
                    old_val = old_row.get(col, "")
                    new_val = new_row.get(col, "")
                    if old_val != new_val:
                        hgvs_diff_rows.extend(
                            component_diffs(urn, col, old_val, new_val)
                        )

    # Write (a) changed rows
    with open(changed_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=new_headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(changed_rows)
    print(f"Wrote {len(changed_rows)} changed rows → {changed_path}", file=sys.stderr)

    # Write (b) HGVS component diffs
    diff_headers = ["variant_urn", "column", "index", "old_value", "new_value"]
    with open(diffs_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=diff_headers, delimiter="\t")
        writer.writeheader()
        writer.writerows(hgvs_diff_rows)
    print(
        f"Wrote {len(hgvs_diff_rows)} HGVS component diffs → {diffs_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
