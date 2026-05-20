# flatten_dna_variants

`src/flatten_dna_variants.py` — Step 12 of the variant-annotation pipeline (optional).

Expands pipe-delimited DNA candidates to one row per candidate. Earlier pipeline steps (particularly `reverse_translate_protein_variants`) can produce rows where a single protein variant maps to multiple genomic DNA candidates, stored as pipe-delimited values across several columns. This script separates those into individual rows.

Protein-only rows — rows where all DNA variant columns are empty (i.e. variants with no reverse translation) — are dropped from the output.

---

## Behaviour

### What gets expanded

Columns are expanded (split on `|`) when they are:

- **Hard-coded DNA variant columns:** `mapped_hgvs_g`, `mapped_hgvs_c`, `dna_clingen_allele_id`, `mapped_hgvs_g_chromosome`, `mapped_hgvs_g_start`, `mapped_hgvs_g_stop`, `mapped_hgvs_g_ref`, `mapped_hgvs_g_alt`, `mapped_hgvs_c_transcript`, `mapped_hgvs_c_start`, `mapped_hgvs_c_stop`, `mapped_hgvs_c_ref`, `mapped_hgvs_c_alt`, `touches_intronic_region`, `spans_intron`, `reverse_translation_warnings`
- **Any column whose name begins with an annotation prefix:** `alphamissense.`, `clingen_evidence_repository.`, `clinvar.`, `gnomad.`, `mutpred2.`, `revel.`, `spliceai.`

These prefixes cover all annotation columns written by earlier pipeline steps. Columns from `annotate_predictors` (`revel.*`, `alphamissense.*`, `mutpred2.*`) are pipe-aligned to candidates by that step and are correctly expanded here.

### What gets repeated (not expanded)

All other columns — e.g. `raw_hgvs_pro`, `gene_symbol`, `score`, `variant_urn`, VEP columns (`vep.*`), MaveDB columns (`mavedb.*`), and assay metadata — are copied as-is to each output row. For `annotate_predictors`, note that `mutpred2.score` is a single value (not pipe-delimited) even for multi-candidate rows, so it is simply repeated.

### Row filtering

A row is kept if at least one DNA variant column is non-empty. Rows where all DNA columns are empty (protein-only variants without a reverse translation) are silently dropped.

---

## Usage

```bash
# Auto-detect DNA variant columns (recommended)
src/scripts/run_flatten_dna_variants.sh annotated.tsv flattened.tsv

# Explicit column list (override auto-detection)
src/scripts/run_flatten_dna_variants.sh annotated.tsv flattened.tsv \
  --dna-variant-columns mapped_hgvs_g,mapped_hgvs_c,dna_clingen_allele_id
```

---

## CLI options

| Option | Default | Description |
|---|---|---|
| `--dna-variant-columns COL,...` | auto-detected | Comma-separated list of columns to expand on `\|`; overrides auto-detection |
| `--skip N` | `0` | Skip first N data rows |
| `--limit N` | no limit | Process at most N rows (after skip) |

Input and output are always TSV. There is no delimiter option.

---

## Dependencies

- **`pandas`** — included in project dependencies.

---

## Troubleshooting

**`No rows with DNA variants found`**

All rows in the input were protein-only (no reverse translation). Confirm that `reverse_translate_protein_variants` ran successfully and produced `mapped_hgvs_g` / `mapped_hgvs_c` values.

**Annotation columns not being expanded**

Auto-detection covers any column beginning with the known annotation prefixes. If you have custom annotation columns that use different names, supply them explicitly with `--dna-variant-columns`.

**Output has more rows than expected**

Each pipe-delimited candidate in any expanded column becomes a row. A row with 3 DNA candidates produces 3 output rows. If a column has fewer values than the maximum number of candidates in a row, the extra positions are filled with empty strings.
