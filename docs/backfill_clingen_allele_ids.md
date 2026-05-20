# backfill_clingen_allele_ids

## Purpose

Re-queries the ClinGen Allele Registry for rows where the `clingen_allele_id`
column is blank after Step 1 (`map_variants`).

`map_variants` populates `clingen_allele_id` with a protein-level ClinGen ID
for each variant. Occasionally some rows are left blank because of transient
API errors, rate-limit pauses, or network timeouts during the initial run.
This utility fills those gaps without re-running the entire mapping step.

**This is not a substitute for Step 3** (`add_dna_clingen_allele_ids`).  Step 3
resolves *DNA-level* allele IDs for every reverse-translated candidate and
writes them to the separate `dna_clingen_allele_id` column. This script only
concerns the single-valued, protein-level `clingen_allele_id` written by
Step 1.

Apply this script to **Step-1 output only** — before reverse translation adds
pipe-delimited candidate lists to the HGVS columns. If the file has already
been through Step 2, the pipe-delimited HGVS values would be submitted as a
single query string, which will fail.

## Input and output columns

**Input:** Any TSV/CSV produced by `map_variants`, with `clingen_allele_id`
partially populated (rows where it is already filled are left untouched).

| Column used | Description |
|---|---|
| `mapped_hgvs_c` | Transcript HGVS (c./n. format) — queried first |
| `mapped_hgvs_g` | Genomic HGVS (g. format) — fallback |
| `mapped_hgvs_p` | Protein HGVS (p. format) — last resort |
| `clingen_allele_id` | Existing value; rows where this is non-empty are skipped |

**Output:** Same columns as input, with blank `clingen_allele_id` cells
filled where ClinGen returned a valid allele record.

## HGVS priority

For each blank row, the script uses the first non-empty value from:

1. `mapped_hgvs_c` (transcript-level, preferred)
2. `mapped_hgvs_g` (genomic)
3. `mapped_hgvs_p` (protein — least likely to resolve a ClinGen record)

This order matches the priority used by `map_variants` itself, so the
resulting `clingen_allele_id` values are consistent.

ClinGen responses that contain blank-node-style placeholder identifiers
(`_:PA...` / `_:CA...`) are treated as misses and left blank.

## Command

```bash
src/scripts/run_backfill_clingen_allele_id.sh input.tsv output.tsv
```

With concurrency tuning (reduce if hitting ClinGen rate limits):

```bash
src/scripts/run_backfill_clingen_allele_id.sh input.tsv output.tsv \
    --concurrency 3
```

Process only a slice of the file (useful when resuming a failed run):

```bash
src/scripts/run_backfill_clingen_allele_id.sh input.tsv output.tsv \
    --skip 2000 --limit 1000
```

## Options

| Flag | Default | Description |
|---|---|---|
| `--clingen-allele-id COLNAME` | `clingen_allele_id` | Column to populate |
| `--hgvs-c-col COLNAME` | `mapped_hgvs_c` | Source column for c./n. HGVS |
| `--hgvs-g-col COLNAME` | `mapped_hgvs_g` | Source column for g. HGVS |
| `--hgvs-p-col COLNAME` | `mapped_hgvs_p` | Source column for p. HGVS |
| `--concurrency N` | `5` | Max concurrent ClinGen API requests |
| `--max-retries N` | `3` | Retries per request (exponential back-off) |
| `--skip N` | `0` | Skip first N data rows |
| `--limit N` | _(all)_ | Process at most N rows after `--skip` |
| `--log-level LEVEL` | `INFO` | Logging verbosity |

## Notes

- Rows with an existing `clingen_allele_id` value are never overwritten, so the
  script is safe to re-run.
- Deduplicates HGVS strings before querying: if multiple blank rows share the
  same HGVS, only one ClinGen request is made.
- Rate-limit responses (HTTP 429) trigger an exponential back-off retry.
- The script runs inside the `map-variants` Docker container using the same
  image and environment as Step 1. Docker must be running.
- File paths are interpreted relative to the `/work` container mount
  (backed by `./data` on the host, or `VARIANT_DATA_DIR` if set).
