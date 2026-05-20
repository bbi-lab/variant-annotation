# utilities

`src/utilities.py` — general-purpose TSV manipulation commands, accessed via `src/scripts/run_utilities.sh`.

All commands auto-detect the delimiter from the file extension (`.tsv`/`.txt` → tab; otherwise comma). Input and output delimiters match independently.

---

## Commands

### `filter-columns`

Keep or omit specific columns.

```bash
# Keep only the listed columns (in input order)
run_utilities.sh filter-columns input.tsv output.tsv --keep-col id --keep-col gene

# Omit specific columns, keep everything else
run_utilities.sh filter-columns input.tsv output.tsv --omit-col internal_notes
```

| Option | Description |
|---|---|
| `--keep-col COLUMN` | Column to keep. Repeatable; also accepts comma-separated values. Mutually exclusive with `--omit-col`. |
| `--omit-col COLUMN` | Column to omit. Repeatable; also accepts comma-separated values. Mutually exclusive with `--keep-col`. |
| `--csv-field-size-limit BYTES` | Increase for large fields. |

---

### `rename-columns`

Keep columns with optional renaming, or omit columns by name. Supports `--reorder` to output columns in the order given.

```bash
# Keep and rename
run_utilities.sh rename-columns input.tsv output.tsv \
  --keep-col old_name:new_name --keep-col gene --reorder

# Omit
run_utilities.sh rename-columns input.tsv output.tsv --omit-col internal_id
```

| Option | Description |
|---|---|
| `--keep-col SRC[:DEST]` | Column to keep, optionally renamed. Repeatable/comma-separated. |
| `--omit-col COLUMN` | Column to omit. Mutually exclusive with `--keep-col`. |
| `--reorder` | Output columns in the order given by `--keep-col`. Only applies in keep mode. |
| `--csv-field-size-limit BYTES` | Increase for large fields. |

---

### `reorder-columns`

Output columns in a specified order; columns not listed are appended at the end.

```bash
run_utilities.sh reorder-columns input.tsv output.tsv --column-order id,gene,score
```

| Option | Description |
|---|---|
| `--column-order COLUMN,...` | Desired column order. Repeatable/comma-separated. |

---

### `filter-rows`

Keep rows where specified columns have non-empty (or blank) values.

```bash
# Keep rows where any of the listed columns is non-empty
run_utilities.sh filter-rows input.tsv output.tsv --value-col mapping_error --match any

# Keep rows where mapped_hgvs_g, mapped_hgvs_c, or mapped_hgvs_p is non-empty
run_utilities.sh filter-rows input.tsv output.tsv \
  --value-col "mapped_hgvs_g,mapped_hgvs_c,mapped_hgvs_p" --match any

# Keep rows where all listed columns are blank (invert)
run_utilities.sh filter-rows input.tsv output.tsv \
  --value-col mapping_error --value-state blank
```

| Option | Default | Description |
|---|---|---|
| `--value-col COLUMN` | (required) | Column to test. Repeatable/comma-separated. |
| `--match any\|all` | `any` | Require any or all columns to satisfy the condition. |
| `--value-state non-empty\|blank` | `non-empty` | Whether the column should be non-empty or blank. |

---

### `replace-rows`

Merge multiple files, later files overriding rows from earlier files on a composite key. New keys from later files are appended.

```bash
run_utilities.sh replace-rows output.tsv base.tsv patch.tsv --key-col variant_urn
```

| Option | Description |
|---|---|
| `--key-col COLUMN` | Key column(s) for row identity. Repeatable/comma-separated. |

---

### `merge-columns`

Left-join a base file with an extra file and add selected columns from the extra file. Supports renaming added columns. The extra file may be a TSV or an Excel workbook (`.xlsx`/`.xls`).

```bash
# Add a single column
run_utilities.sh merge-columns base.tsv extra.tsv output.tsv \
  --key-col variant_urn --add-col score

# Add and rename; key columns have different names in the two files
run_utilities.sh merge-columns base.tsv extra.tsv output.tsv \
  --key-col "dataset_name:Dataset Name" \
  --add-col "Gene:gene_symbol" \
  --add-col "HGNC ID:hgnc_id"

# Read a specific worksheet from an Excel extra file
run_utilities.sh merge-columns base.tsv metadata.xlsx output.tsv \
  --key-col id --add-col notes --extra-worksheet Curation
```

| Option | Description |
|---|---|
| `--key-col BASE_COL[:EXTRA_COL]` | Join key. Use `BASE:EXTRA` when columns have different names. Repeatable/comma-separated. |
| `--add-col SRC[:DEST]` | Column to add from extra file, optionally renamed. Repeatable/comma-separated. |
| `--add-all-cols-from-extra` | Add all non-key columns from extra not already in base. |
| `--extra-worksheet NAME\|INDEX` | Worksheet to read from an Excel extra file (name or 1-based index). |

---

### `compare-columns`

Find rows where paired column values differ. Writes differing rows with a `differences` column listing the mismatched pairs. Output goes to `--output` or stdout.

```bash
run_utilities.sh compare-columns input.tsv \
  --col-a mapped_hgvs_g --col-b ref_hgvs_g \
  --output diffs.tsv
```

| Option | Description |
|---|---|
| `--col-a COLUMN` | Column for group A. Repeatable to compare multiple pairs. |
| `--col-b COLUMN` | Column for group B. Must match the count of `--col-a`. |
| `--output FILE` | Write output to FILE (default: stdout). |
| `--skip N` | Skip first N rows. |
| `--limit N` | Process at most N rows. |
