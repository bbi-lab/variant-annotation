"""
Compare the VCF-style coordinate columns between two TSV files, joining on variant_urn.

For genome, transcript, and protein groups, compares the pipe-separated components of
chromosome/transcript, start, stop, ref, and alt.

Groups:
  genome     — mapped_hgvs_g and mapped_hgvs_g_{chromosome,start,stop,ref,alt}
  transcript — mapped_hgvs_c and mapped_hgvs_c_{transcript,start,stop,ref,alt}
  protein    — mapped_hgvs_p and mapped_hgvs_p_{chromosome,start,stop,ref,alt}

A difference in any sub-column within the same group and index position yields one row:
  variant_urn, group, index,
  old_hgvs, new_hgvs,
  old_chromosome (or old_transcript for the transcript group), new_chromosome/new_transcript,
  old_start, new_start, old_stop, new_stop, old_ref, new_ref, old_alt, new_alt

If a variant_urn exists in only one file, the missing side's values are blank.
Groups whose columns are absent from both files are silently skipped.

Output defaults to <new_tsv_stem>.vcf_diffs.tsv in the current directory.
Override with -o / --output.

Both files must contain a variant_urn column.
Runs locally using the project virtual environment (not Docker).
"""

import argparse
import csv
import sys
from pathlib import Path

csv.field_size_limit(10 * 1024 * 1024)  # 10 MB

KEY_COL = "variant_urn"

# Each group: (group_label, hgvs_col, locus_col, locus_header, start_col, stop_col, ref_col, alt_col)
VCF_GROUPS = [
    (
        "genome",
        "mapped_hgvs_g",
        "mapped_hgvs_g_chromosome",
        "chromosome",
        "mapped_hgvs_g_start",
        "mapped_hgvs_g_stop",
        "mapped_hgvs_g_ref",
        "mapped_hgvs_g_alt",
    ),
    (
        "transcript",
        "mapped_hgvs_c",
        "mapped_hgvs_c_transcript",
        "transcript",
        "mapped_hgvs_c_start",
        "mapped_hgvs_c_stop",
        "mapped_hgvs_c_ref",
        "mapped_hgvs_c_alt",
    ),
    (
        "protein",
        "mapped_hgvs_p",
        "mapped_hgvs_p_chromosome",
        "chromosome",
        "mapped_hgvs_p_start",
        "mapped_hgvs_p_stop",
        "mapped_hgvs_p_ref",
        "mapped_hgvs_p_alt",
    ),
]

ALL_VCF_COLS = {col for group in VCF_GROUPS for col in (group[1],) + group[2:3] + group[4:]}
HGVS_COLS = {"mapped_hgvs_g", "mapped_hgvs_c", "mapped_hgvs_p"}

DIFF_HEADERS = [
    "variant_urn",
    "group",
    "index",
    "old_hgvs",
    "new_hgvs",
    "old_chromosome",
    "new_chromosome",
    "old_start",
    "new_start",
    "old_stop",
    "new_stop",
    "old_ref",
    "new_ref",
    "old_alt",
    "new_alt",
]


def load_cols(path: str) -> dict[str, dict[str, str]]:
    """Return {variant_urn: {col: value}} for all VCF coordinate and HGVS columns present."""
    rows: dict[str, dict[str, str]] = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise ValueError(f"Empty file: {path}")
        present = (set(reader.fieldnames) & ALL_VCF_COLS) | (set(reader.fieldnames) & HGVS_COLS)
        for row in reader:
            key = row[KEY_COL]
            rows[key] = {col: row.get(col, "") for col in present}
    return rows


def group_diffs(
    urn: str,
    group_label: str,
    hgvs_col: str,
    old_row: dict[str, str],
    new_row: dict[str, str],
    locus_col: str,
    start_col: str,
    stop_col: str,
    ref_col: str,
    alt_col: str,
) -> list[dict[str, str]]:
    """Compare pipe-separated components for a group; yield one row per changed index."""

    def parts(row: dict[str, str], col: str) -> list[str]:
        val = row.get(col, "")
        return val.split("|") if val else []

    old_hgvs  = parts(old_row, hgvs_col)
    old_locus = parts(old_row, locus_col)
    old_start = parts(old_row, start_col)
    old_stop  = parts(old_row, stop_col)
    old_ref   = parts(old_row, ref_col)
    old_alt   = parts(old_row, alt_col)

    new_hgvs  = parts(new_row, hgvs_col)
    new_locus = parts(new_row, locus_col)
    new_start = parts(new_row, start_col)
    new_stop  = parts(new_row, stop_col)
    new_ref   = parts(new_row, ref_col)
    new_alt   = parts(new_row, alt_col)

    length = max(
        len(old_locus), len(old_start), len(old_stop), len(old_ref), len(old_alt),
        len(new_locus), len(new_start), len(new_stop), len(new_ref), len(new_alt),
    )
    if length == 0:
        return []

    def get(lst: list[str], i: int) -> str:
        return lst[i] if i < len(lst) else ""

    diffs = []
    for i in range(length):
        ol, nl = get(old_locus, i), get(new_locus, i)
        os_, ns = get(old_start, i), get(new_start, i)
        oe, ne = get(old_stop,  i), get(new_stop,  i)
        or_, nr = get(old_ref,  i), get(new_ref,   i)
        oa, na = get(old_alt,   i), get(new_alt,   i)
        if ol != nl or os_ != ns or oe != ne or or_ != nr or oa != na:
            diffs.append(
                {
                    "variant_urn": urn,
                    "group": group_label,
                    "index": str(i),
                    "old_hgvs":       get(old_hgvs, i),
                    "new_hgvs":       get(new_hgvs, i),
                    "old_chromosome": ol,
                    "new_chromosome": nl,
                    "old_start": os_,
                    "new_start": ns,
                    "old_stop":  oe,
                    "new_stop":  ne,
                    "old_ref":   or_,
                    "new_ref":   nr,
                    "old_alt":   oa,
                    "new_alt":   na,
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
        help="Output file path (default: <new_tsv_stem>.vcf_diffs.tsv)",
    )
    args = parser.parse_args()

    output_path = args.output or f"{Path(args.new_tsv).stem}.vcf_diffs.tsv"

    print(f"Loading {args.old_tsv} ...", file=sys.stderr)
    old_rows = load_cols(args.old_tsv)
    print(f"Loading {args.new_tsv} ...", file=sys.stderr)
    new_rows = load_cols(args.new_tsv)

    all_urns = old_rows.keys() | new_rows.keys()
    empty: dict[str, str] = {}

    diff_rows: list[dict[str, str]] = []
    for urn in sorted(all_urns):
        old = old_rows.get(urn, empty)
        new = new_rows.get(urn, empty)
        for group_label, hgvs_col, locus_col, _locus_header, start_col, stop_col, ref_col, alt_col in VCF_GROUPS:
            diff_rows.extend(
                group_diffs(
                    urn, group_label, hgvs_col, old, new,
                    locus_col, start_col, stop_col, ref_col, alt_col,
                )
            )

    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=DIFF_HEADERS, delimiter="\t")
        writer.writeheader()
        writer.writerows(diff_rows)

    print(f"Wrote {len(diff_rows)} VCF coordinate diffs → {output_path}", file=sys.stderr)


if __name__ == "__main__":
    main()
