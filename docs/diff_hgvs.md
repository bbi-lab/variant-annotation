# diff_hgvs

`src/diff_hgvs.py` — compare the HGVS mapping columns between two versions of a variants TSV file.

Joins old and new files on `variant_urn` and reports every component of `mapped_hgvs_g`, `mapped_hgvs_c`, and `mapped_hgvs_p` that changed. Each changed pipe-delimited component produces one output row.

Runs locally using the project's Python virtual environment (not Docker).

---

## Output

One TSV row per changed component:

| Column | Description |
|---|---|
| `variant_urn` | Variant identifier |
| `column` | Which HGVS column changed (`mapped_hgvs_g`, `mapped_hgvs_c`, or `mapped_hgvs_p`) |
| `index` | 0-based position within the pipe-delimited value |
| `old_value` | Component value in the old file (blank if not present) |
| `new_value` | Component value in the new file (blank if not present) |

Variants present in only one file appear with the missing side's components left blank.

Output defaults to `<new_tsv_stem>.hgvs_diffs.tsv` in the current directory.

---

## Usage

```bash
src/scripts/run_diff_hgvs.sh old.tsv new.tsv
src/scripts/run_diff_hgvs.sh old.tsv new.tsv output.tsv
```

| Argument | Description |
|---|---|
| `old_tsv` | Baseline TSV file |
| `new_tsv` | Updated TSV file |
| `output.tsv` | Optional output path (third positional argument) |
