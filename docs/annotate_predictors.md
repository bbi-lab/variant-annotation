# annotate_predictors

`src/annotate_predictors.py` — Step 11 of the variant-annotation pipeline (optional).

Annotates each variant with pre-computed in-silico pathogenicity scores from computational predictors. All currently supported predictors score **missense SNVs only** (GRCh38); other variant types (indels, synonymous substitutions, splice-site variants, etc.) receive empty annotation columns.

At least one of `--revel-file`, `--alphamissense-file`, or `--dbnsfp-file` must be supplied.

---

## Supported tools and output columns

| Tool | Output column(s) | Range | Column flag | Env variable |
|---|---|---|---|---|
| REVEL | `revel.score` | 0–1 (higher = more pathogenic) | `--revel-file` | `REVEL_FILE` |
| AlphaMissense | `alphamissense.pathogenicity` | 0–1 | `--alphamissense-file` | `ALPHAMISSENSE_FILE` |
| AlphaMissense | `alphamissense.class` | `likely_benign` / `ambiguous` / `likely_pathogenic` | (same file) | |
| MutPred2 (via dbNSFP) | `mutpred2.score` | 0–1 | `--dbnsfp-file` | `DBNSFP_FILE` |

Only the columns for tools whose data files are provided will appear in the output. All score values are formatted to 4 decimal places.

---

## Pipe-delimited column behaviour

When a row has multiple DNA candidates (pipe-delimited `mapped_hgvs_g`):

- **REVEL** and **AlphaMissense** columns are pipe-aligned to the input candidates — each candidate gets its own score value (or empty if not an SNV or not found).
- **MutPred2** is a protein-level model. All reverse-translation candidates encode the same amino acid substitution, so a single maximum score across all SNV candidates is emitted (not pipe-delimited).

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

### MutPred2 (via dbNSFP)

1. Download the pre-built GRCh38 variant file and its `.tbi` index from:
   `https://sites.google.com/site/jpopgen/dbNSFP`
   (e.g. `dbNSFP5.3.1a_grch38.gz` and `dbNSFP5.3.1a_grch38.gz.tbi`)
2. Both files must be in the same directory. No further preparation required.

---

## Usage

```bash
# AlphaMissense only
src/scripts/run_annotate_predictors.sh input.tsv output.tsv \
  --alphamissense-file AlphaMissense_hg38.tsv.gz

# All three tools
src/scripts/run_annotate_predictors.sh input.tsv output.tsv \
  --revel-file revel_hg38.tsv.gz \
  --alphamissense-file AlphaMissense_hg38.tsv.gz \
  --dbnsfp-file dbNSFP5.3.1a_grch38.gz
```

---

## CLI options

| Option | Env variable | Default | Description |
|---|---|---|---|
| `--revel-file PATH` | `REVEL_FILE` | — | bgzipped, tabix-indexed REVEL TSV |
| `--alphamissense-file PATH` | `ALPHAMISSENSE_FILE` | — | bgzipped, tabix-indexed AlphaMissense TSV |
| `--dbnsfp-file PATH` | `DBNSFP_FILE` | — | bgzipped, tabix-indexed dbNSFP GRCh38 variant file |
| `--mapped-hgvs-g-col COL` | — | `mapped_hgvs_g` | Input column containing genomic HGVS values |
| `--skip N` | — | `0` | Skip first N data rows |
| `--limit N` | — | no limit | Stop after N rows |
| `--log-level` | — | `INFO` | Logging verbosity |
| `--csv-field-size-limit BYTES` | — | system default | Increase for large fields |

The input delimiter is auto-detected from the file extension (`.tsv`/`.txt` → tab; otherwise comma). Output delimiter matches.

---

## Dependencies

- **`tabix`** (htslib) — required system-level tool. The script verifies `tabix` is available at startup and exits with an error if not found.
- No additional Python library dependencies.

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

dbNSFP stores multiple MutPred2 scores per row (one per transcript), semicolon-separated. The script returns the maximum non-null value across all entries for each allele.
