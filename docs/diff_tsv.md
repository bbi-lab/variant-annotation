# diff_tsv

`src/diff_tsv.py` — row-level diff between two versions of a variants TSV file.

Joins old and new files on `variant_urn` and produces two output files: one with all rows from the new file that differ in any column from the old file, and one with per-component HGVS diffs (equivalent to `diff_hgvs`).

Runs locally using the project's Python virtual environment (not Docker).

---

## Output files

| File | Content |
|---|---|
| `<prefix>.changed_rows.tsv` | All rows from the new file whose content differs from the corresponding row in the old file; includes only rows present in both files whose values changed. |
| `<prefix>.hgvs_diffs.tsv` | Per-component diffs for `mapped_hgvs_g`, `mapped_hgvs_c`, and `mapped_hgvs_p` (same format as `diff_hgvs` output: `variant_urn, column, index, old_value, new_value`). |

The output prefix defaults to the stem of the new TSV filename.

---

## Usage

```bash
src/scripts/run_diff_tsv.sh old.tsv new.tsv
src/scripts/run_diff_tsv.sh old.tsv new.tsv output_prefix
```

| Argument | Description |
|---|---|
| `old_tsv` | Baseline TSV file |
| `new_tsv` | Updated TSV file |
| `output_prefix` | Optional prefix for output file names (third positional argument) |

---

## Relationship to other diff tools

`diff_tsv` is a combined tool: its `.changed_rows.tsv` output shows row-level changes across all columns, while its `.hgvs_diffs.tsv` output is equivalent to running `diff_hgvs` separately. Use `diff_hgvs` or `diff_vcf_coords` directly when you only need coordinate-specific diffs.
