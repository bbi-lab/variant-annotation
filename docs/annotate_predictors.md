# annotate_predictors

`src/annotate_predictors.py` — Step 11 of the variant-annotation pipeline (optional).

Annotates each variant with pre-computed in-silico pathogenicity scores from computational predictors. All currently supported predictors score **missense SNVs only** (GRCh38); other variant types (indels, synonymous substitutions, splice-site variants, etc.) receive empty annotation columns.

Tabix queries (REVEL, AlphaMissense, dbNSFP) use thread-local `pysam.TabixFile` handles when pysam is available, falling back to subprocess `tabix` otherwise. `--max-workers N` (N > 1) parallelizes those lookups across threads within a single process: a first pass collects all unique SNVs in the input, looks them up concurrently, then a second pass streams the annotated output.

At least one of `--revel-file`, `--revel-cache-file`, `--alphamissense-file`, `--alphamissense-cache-file`, `--mutpred2-properties-file`, `--dbnsfp-file`, `--revel-training-file`, or `--mutpred2-training-file` must be supplied.

---

## Supported tools and output columns

| Tool | Output column(s) | Range | Column flag | Env variable |
|---|---|---|---|---|
| REVEL | `revel.score` | 0–1 (higher = more pathogenic) | `--revel-file` | `REVEL_FILE` |
| AlphaMissense | `alphamissense.pathogenicity` | 0–1 | `--alphamissense-file` | `ALPHAMISSENSE_FILE` |
| AlphaMissense | `alphamissense.class` | `likely_benign` / `ambiguous` / `likely_pathogenic` | (same file) | |
| MutPred2 (via properties file — preferred) | `mutpred2.score` | 0–1 | `--mutpred2-properties-file` | `MUTPRED2_PROPERTIES_FILE` |
| MutPred2 (via dbNSFP — legacy fallback) | `mutpred2.score` | 0–1 | `--dbnsfp-file` | `DBNSFP_FILE` |
| REVEL training-set overlap | `revel.train` | `true` / `false` | `--revel-training-file` | `REVEL_TRAINING_FILE` |
| MutPred2 training-set overlap | `mutpred2.train` | `true` / `false` | `--mutpred2-training-file` | `MUTPRED2_TRAINING_FILE` |

Only the columns for tools whose data files are provided will appear in the output. All score values are formatted to 4 decimal places.

MutPred2 has two, mutually exclusive sources. If `--mutpred2-properties-file` is given it takes precedence and `--dbnsfp-file` is ignored (a warning is logged when both are supplied).

The two training-set columns are independent, optional additions — they can be supplied with or without any of the score sources above, and don't require `tabix`.

---

## Pipe-delimited column behaviour

When a row has multiple DNA candidates (pipe-delimited `mapped_hgvs_g`):

- **REVEL** and **AlphaMissense** columns are pipe-aligned to the input candidates — each candidate gets its own score value (or empty if not an SNV or not found).
- **MutPred2 (properties file)** is looked up per reverse-translation candidate, keyed on `variant_urn` plus the genomic coordinates for that candidate (chromosome, start, stop, ref, alt), so it is pipe-aligned like REVEL/AlphaMissense. A candidate absent from the properties file produces an empty slot rather than being dropped.
- **MutPred2 (dbNSFP, legacy)** is treated as a protein-level model there: all reverse-translation candidates encode the same amino acid substitution, so a single maximum score across all SNV candidates is emitted (not pipe-delimited).
- **`revel.train`** is looked up per candidate by genomic coordinates, like REVEL/AlphaMissense scores, and is `true`/`false` per pipe slot.
- **`mutpred2.train`** is protein-level like MutPred2 via dbNSFP, but for consistency with the other pipe-aligned columns it is still emitted once per DNA candidate — the same `true`/`false` value is duplicated across every pipe slot in the row.

---

## File-based caches

For large batch jobs or offline use, REVEL and AlphaMissense scores can be pre-loaded from TSV files rather than queried via tabix. File-based caches are checked first; tabix is used as a fallback for any variant not found in the cache. When only file caches are configured (no tabix files), `tabix` does not need to be installed.

### REVEL cache file

A two-column tab-separated file (header required):

| Column | Description |
|---|---|
| `hgvs` | Transcript HGVS c-string (e.g. `NM_000133.4:c.1A>G`) |
| `revel.score` | Pre-computed REVEL score |

Passed via `--revel-cache-file PATH` or the `REVEL_CACHE_FILE` environment variable.

### AlphaMissense cache file

A three-column tab-separated file (header required):

| Column | Description |
|---|---|
| `hgvs` | Transcript HGVS c-string |
| `alphamissense.pathogenicity` | Pre-computed pathogenicity score |
| `alphamissense.class` | Pre-computed class label |

Passed via `--alphamissense-cache-file PATH` or the `ALPHAMISSENSE_CACHE_FILE` environment variable.

### Cache keys and fallback

Both caches are keyed on transcript HGVS c-strings from the column named by `--mapped-hgvs-c-col` (default: `mapped_hgvs_c`). If the c-string for a candidate is blank, the lookup falls back to tabix using the genomic HGVS from `--mapped-hgvs-g-col`. File-cache keys and tabix lookups are therefore complementary: a cache hit skips the tabix query for that candidate, while uncached variants still hit tabix normally.

A typical use case is to extract previously computed scores from an existing output TSV, de-duplicate them by HGVS key, and feed the result back in as a warm cache to skip redundant tabix queries on incremental runs.

---

## Training-set overlap files

These flag whether each variant was part of the training set for a predictor, so held-out performance can be assessed separately from variants the model was trained on. Both are optional, independent of the score sources above, and don't require `tabix`.

### REVEL training-variants file

A tab-separated file (header required) with at least these columns:

| Column | Description |
|---|---|
| `chromosome` | Chromosome (e.g. `17`) |
| `hg38_start` | Genomic start position (hg38); floats like `3476169.0` are tolerated |
| `hg38_end` | Genomic end position (hg38) |
| `ref_allele` | Reference allele |
| `alt_allele` | Alternate allele |

Any other columns (e.g. `gene_symbol`, `unqualified_hgvs_p`) are informational and ignored for matching. See `data/revel_training_variants.tsv` for an example.

Passed via `--revel-training-file PATH` or the `REVEL_TRAINING_FILE` environment variable. Joined against the `--mapped-hgvs-g-chromosome-col`/`-start-col`/`-stop-col`/`-ref-col`/`-alt-col` columns (the same columns used for `--mutpred2-properties-file`), per DNA reverse-translation candidate. Produces `revel.train`, pipe-aligned like `revel.score`. Chromosome is tried with and without a `chr` prefix; positions and alleles must match verbatim.

### MutPred2 training-variants file

A tab-separated file (header required) using one of two schemas, auto-detected from the header:

| Schema | Required column(s) | Notes |
|---|---|---|
| Qualified | `hgvs_p` | Transcript/protein-accession-qualified, e.g. `NP_000546.1:p.Asn1643His`. Matched verbatim against `--mapped-hgvs-p-col`; gene symbols are ignored. |
| Unqualified | `gene_symbol`, `unqualified_hgvs_p` | Matched against `--gene-symbol-col` and the part of `--mapped-hgvs-p-col` following the colon (e.g. `p.Asn1643His`). Used when there is no `hgvs_p` column. |

See `data/mutpred2_training_variants.tsv` for an example (unqualified schema).

Passed via `--mutpred2-training-file PATH` or the `MUTPRED2_TRAINING_FILE` environment variable. Because MutPred2 is a protein-level model, the match is evaluated once per row and the same `true`/`false` value is then duplicated across every pipe slot in `mutpred2.train`, so it still aligns with the other DNA-candidate columns.

If `--mapped-hgvs-p-col` itself holds a pipe-delimited value (e.g. multiple transcript mappings for one protein change), a match on any segment counts as a match for the whole row.

---

## Data file preparation

### AlphaMissense

1. Download `AlphaMissense_hg38.tsv.gz` from:
   `https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz`
2. Generate the tabix index:
   ```bash
   tabix -s 1 -b 2 -e 2 AlphaMissense_hg38.tsv.gz
   ```

### REVEL

1. Download `revel_with_transcript_ids.csv.zip` from:
   `https://sites.google.com/site/revelgenomics/downloads`
2. Prepare the indexed file:
   ```bash
   unzip revel_with_transcript_ids.csv.zip
   tail -n +2 revel_with_transcript_ids.csv \
     | awk -F',' 'NF>=9 && $3!="" && $3!="." \
                  {print $1"\t"$3"\t"$4"\t"$5"\t"$8}' \
     | (printf '#chr\tpos\tref\talt\trevel_score\n'; sort -k1,1V -k2,2n) \
     | bgzip > revel_hg38.tsv.gz
   tabix -s 1 -b 2 -e 2 -S 1 revel_hg38.tsv.gz
   ```
   The resulting file has five tab-separated columns: `#chr pos ref alt revel_score`.

### MutPred2 (via properties file)

A MaveDB-derived per-DNA-variant properties table (CSV, optionally gzipped) with (at least) these columns:

| Column | Description |
|---|---|
| `mavedb_variant_urn` | MaveDB variant URN (protein- or DNA-level) |
| `Chrom` | Chromosome (with or without `chr` prefix — both are tried) |
| `hg38_start` | Genomic start position (hg38) |
| `hg38_end` | Genomic end position (hg38) |
| `ref_allele` | Reference allele |
| `alt_allele` | Alternate allele |
| `MutPred2 score` | Pre-computed MutPred2 score |

The same `mavedb_variant_urn` can repeat across multiple rows (one per DNA reverse-translation candidate); each is keyed separately on `(variant_urn, chrom, start, stop, ref, alt)`. No further preparation required beyond having the file on disk.

### MutPred2 (via dbNSFP, legacy)

1. Download the pre-built GRCh38 variant file and its `.tbi` index from:
   `https://sites.google.com/site/jpopgen/dbNSFP`
   (e.g. `dbNSFP5.3.1a_grch38.gz` and `dbNSFP5.3.1a_grch38.gz.tbi`)
2. Both files must be in the same directory. No further preparation required.
3. Ignored if `--mutpred2-properties-file` is also given.

---

## Usage

```bash
# AlphaMissense only
src/scripts/run_annotate_predictors.sh input.tsv output.tsv \
  --alphamissense-file AlphaMissense_hg38.tsv.gz

# All three tools, MutPred2 via legacy dbNSFP source
src/scripts/run_annotate_predictors.sh input.tsv output.tsv \
  --revel-file revel_hg38.tsv.gz \
  --alphamissense-file AlphaMissense_hg38.tsv.gz \
  --dbnsfp-file dbNSFP5.3.1a_grch38.gz

# All three tools, MutPred2 via preferred properties-file source (pipe-aligned)
src/scripts/run_annotate_predictors.sh input.tsv output.tsv \
  --revel-file revel_hg38.tsv.gz \
  --alphamissense-file AlphaMissense_hg38.tsv.gz \
  --mutpred2-properties-file data_frame_missense_variants_MP2_properties.csv.gz

# Training-set overlap flags only, no scores
src/scripts/run_annotate_predictors.sh input.tsv output.tsv \
  --revel-training-file data/revel_training_variants.tsv \
  --mutpred2-training-file data/mutpred2_training_variants.tsv
```

---

## CLI options

| Option | Env variable | Default | Description |
|---|---|---|---|
| `--revel-file PATH` | `REVEL_FILE` | — | bgzipped, tabix-indexed REVEL TSV |
| `--revel-cache-file PATH` | `REVEL_CACHE_FILE` | — | Two-column TSV `(hgvs, revel.score)` pre-loaded as a file-based REVEL cache; checked before tabix |
| `--alphamissense-file PATH` | `ALPHAMISSENSE_FILE` | — | bgzipped, tabix-indexed AlphaMissense TSV |
| `--alphamissense-cache-file PATH` | `ALPHAMISSENSE_CACHE_FILE` | — | Three-column TSV `(hgvs, alphamissense.pathogenicity, alphamissense.class)` pre-loaded as a file-based AlphaMissense cache; checked before tabix |
| `--dbnsfp-file PATH` | `DBNSFP_FILE` | — | bgzipped, tabix-indexed dbNSFP GRCh38 variant file. Ignored if `--mutpred2-properties-file` is also given |
| `--mutpred2-properties-file PATH` | `MUTPRED2_PROPERTIES_FILE` | — | MaveDB MP2-properties CSV (optionally gzipped); preferred source for `mutpred2.score`, takes precedence over `--dbnsfp-file` |
| `--revel-training-file PATH` | `REVEL_TRAINING_FILE` | — | TSV of REVEL training variants (genomic coordinates); produces `revel.train` |
| `--mutpred2-training-file PATH` | `MUTPRED2_TRAINING_FILE` | — | TSV of MutPred2 training variants (`hgvs_p`, or `gene_symbol` + `unqualified_hgvs_p`); produces `mutpred2.train` |
| `--variant-urn-col COL` | — | `variant_urn` | Input column with the MaveDB variant URN, used as part of the lookup key for `--mutpred2-properties-file` |
| `--mapped-hgvs-g-chromosome-col COL` | — | `mapped_hgvs_g_chromosome` | Input column with pipe-delimited genomic chromosome(s), used as part of the lookup key for `--mutpred2-properties-file` and `--revel-training-file` |
| `--mapped-hgvs-g-start-col COL` | — | `mapped_hgvs_g_start` | Input column with pipe-delimited genomic start position(s), used as part of the lookup key for `--mutpred2-properties-file` and `--revel-training-file` |
| `--mapped-hgvs-g-stop-col COL` | — | `mapped_hgvs_g_stop` | Input column with pipe-delimited genomic stop position(s), used as part of the lookup key for `--mutpred2-properties-file` and `--revel-training-file` |
| `--mapped-hgvs-g-ref-col COL` | — | `mapped_hgvs_g_ref` | Input column with pipe-delimited genomic ref allele(s), used as part of the lookup key for `--mutpred2-properties-file` and `--revel-training-file` |
| `--mapped-hgvs-g-alt-col COL` | — | `mapped_hgvs_g_alt` | Input column with pipe-delimited genomic alt allele(s), used as part of the lookup key for `--mutpred2-properties-file` and `--revel-training-file` |
| `--gene-symbol-col COL` | — | `gene_symbol` | Input column with the gene symbol, used for `--mutpred2-training-file` when that file has no `hgvs_p` column |
| `--mapped-hgvs-p-col COL` | — | `mapped_hgvs_p` | Input column with protein HGVS value(s), used as part of the lookup key for `--mutpred2-training-file` |
| `--mapped-hgvs-g-col COL` | — | `mapped_hgvs_g` | Input column containing pipe-delimited genomic HGVS values |
| `--mapped-hgvs-c-col COL` | — | `mapped_hgvs_c` | Input column containing pipe-delimited transcript HGVS values used as file-cache keys |
| `--skip N` | — | `0` | Skip first N data rows |
| `--limit N` | — | no limit | Stop after N rows |
| `--log-level` | — | `INFO` | Logging verbosity |
| `--csv-field-size-limit BYTES` | — | system default | Increase for large fields |
| `--max-workers N` | — | `1` | Parallel tabix lookup threads for REVEL/AlphaMissense/dbNSFP. No effect on `--mutpred2-properties-file`, `--revel-training-file`, or `--mutpred2-training-file`, which don't use tabix |

The input delimiter is auto-detected from the file extension (`.tsv`/`.txt` → tab; otherwise comma). Output delimiter matches.

---

## Dependencies

- **`tabix`** (htslib) — required when any tabix-indexed data file is configured (`--revel-file`, `--alphamissense-file`, `--dbnsfp-file`). Not required when running with file-based caches alone. The script verifies `tabix` is available at startup only when needed.
- **`pysam`** (optional but recommended) — enables faster thread-local tabix handles, avoiding subprocess overhead per lookup, and is required for `--max-workers` to actually run lookups in parallel (without it, per-query subprocess spawns serialize under the GIL). If pysam is not installed, every tabix query spawns a subprocess.
- **`pandas`** — used to load `--mutpred2-properties-file` (a multi-column, often multi-million-row CSV) with its C parser rather than a per-row Python loop. Not required for any other code path.

---

## Troubleshooting

**`tabix executable not found`**

Install [htslib](http://www.htslib.org/) and ensure `tabix` is on `$PATH`.

**All score columns are empty**

- Confirm `--mapped-hgvs-g-col` points to a column with genomic HGVS strings in the format `NC_XXXXXX.XX:g.POS REF>ALT`. Non-SNV variants (indels, MNVs) produce empty scores by design.
- The data file covers missense SNVs only. Synonymous, stop-gain, and non-coding variants will always be empty.

**REVEL returns no scores**

- The REVEL preparation step filters for rows where the position field is non-empty (`$3 != ""` and `$3 != "."`). If the file was prepared differently, scores may not align with the expected column layout.

**AlphaMissense: multiple transcripts per allele**

The script returns the highest pathogenicity score when multiple transcript entries exist for the same allele; the associated `am_class` label corresponds to that highest-scoring entry.

**MutPred2: semicolon-separated values**

dbNSFP stores multiple MutPred2 scores per row (one per transcript), semicolon-separated. The script returns the maximum non-null value across all entries for each allele. This only applies to the legacy `--dbnsfp-file` source.

**MutPred2 (properties file): candidate not found**

Confirm `--variant-urn-col` and the `--mapped-hgvs-g-*-col` columns are populated and match the `mavedb_variant_urn` / `Chrom` / `hg38_start` / `hg38_end` / `ref_allele` / `alt_allele` values in the properties file exactly (chromosome is tried with and without a `chr` prefix, but positions and alleles must match verbatim).

**Slow tabix lookups (REVEL/AlphaMissense/dbNSFP)**

Ensure `pysam` is installed so lookups use thread-local `TabixFile` handles instead of spawning a `tabix` subprocess per query — this alone is typically the largest speedup. Then increase `--max-workers` to parallelize those lookups across threads within the process, rather than splitting the input and running multiple separate processes.

**`revel.train` / `mutpred2.train` are always `false`**

- For `revel.train`, confirm the `--mapped-hgvs-g-*-col` columns are populated (they default to the same columns used by `--mutpred2-properties-file`, e.g. produced by `add_vcf_identifiers`) and that positions/alleles match the training file's `hg38_start`/`hg38_end`/`ref_allele`/`alt_allele` verbatim.
- For `mutpred2.train`, check which schema the training file was auto-detected as (logged at startup, e.g. `schema: unqualified`). With the unqualified schema, confirm `--gene-symbol-col` matches the input's gene symbol column and that `--mapped-hgvs-p-col` includes the accession qualifier consistently (e.g. `NP_000546.1:p.Asn1643His`) so the part after the colon lines up with the training file's `unqualified_hgvs_p`.
- The end-of-run log line reports both the row/genomic-variant match counts and which columns were used for the join, which is usually the fastest way to spot a column-name mismatch.
