# annotate_clinvar

`src/annotate_clinvar.py` — Step 5 of the variant-annotation pipeline (optional).

Looks up ClinVar clinical classification data for each variant using the monthly ClinVar variant-summary TSV archive from NCBI. For each row the script resolves the `dna_clingen_allele_id` (CA-prefixed ClinGen identifier) to a ClinVar Allele ID via the ClinGen Allele Registry, then looks up that Allele ID in the cached TSV to retrieve clinical significance, review status, star rating, and last review date.

---

## What it does

1. Downloads (if not already cached) the ClinVar variant-summary TSV for the requested month from the NCBI FTP archive.
2. Parses the TSV into an in-memory dict keyed by ClinVar Allele ID. Only GRCh38 rows are retained.
3. For each input row, resolves the `dna_clingen_allele_id` to a `(ClinVar Variation ID, ClinVar Allele ID)` pair by querying the ClinGen Allele Registry API.
4. Looks up the ClinVar Allele ID in the in-memory dict and writes six annotation columns.
5. Rows are processed concurrently and written in input order.

When `dna_clingen_allele_id` is pipe-delimited (as produced by `add_dna_clingen_allele_ids` for reverse-translated protein rows), each candidate is resolved and annotated independently. Output columns are also pipe-delimited, preserving candidate cardinality.

---

## Output columns

Output column names follow the pattern `<namespace>.<version>.<field>`, where `namespace` defaults to `clinvar` and `version` is the YYYYMM release string.

| Column | Description |
|---|---|
| `clinvar.<YYYYMM>.variation_id` | ClinVar Variation ID (numeric string) for the resolved allele |
| `clinvar.<YYYYMM>.allele_id` | ClinVar Allele ID used to look up the TSV record |
| `clinvar.<YYYYMM>.clinical_significance` | ClinVar `ClinicalSignificance` field (e.g. `Pathogenic`, `Likely benign`) |
| `clinvar.<YYYYMM>.review_status` | ClinVar `ReviewStatus` field (e.g. `criteria provided, single submitter`) |
| `clinvar.<YYYYMM>.stars` | Star rating derived from review status (see [Star ratings](#star-ratings)) |
| `clinvar.<YYYYMM>.last_review_date` | ClinVar `LastEvaluated` date string |

ClinVar `-` placeholder values are normalised to empty strings. When no ClinVar record is found for a candidate, all six fields are empty for that candidate position.

---

## Star ratings

The `stars` column maps the ClinVar review status to a 0–4 integer following ClinVar's own definitions:

| Stars | Review status |
|---|---|
| 4 | practice guideline |
| 3 | reviewed by expert panel |
| 2 | criteria provided, multiple submitters, no conflicts |
| 1 | criteria provided, conflicting classifications |
| 1 | criteria provided, single submitter |
| 0 | no assertion criteria provided |
| 0 | no assertion provided |
| 0 | no classifications from unflagged records |

---

## ClinVar TSV caching

The NCBI variant-summary files are monthly archives and never change after release. The downloaded `.gz` file is stored in `--cache-dir` and reused on all subsequent runs with no TTL. The default cache directory is `$CLINVAR_CACHE_DIR` or `/tmp/clinvar_cache`.

The file for a given `--clinvar-version` is downloaded from:

```
https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/variant_summary_{YYYY}-{MM}.txt.gz
```

with a fallback to the year-subdirectory URL path if the primary URL returns an error.

---

## ClinGen → ClinVar ID resolution

Resolving a ClinGen allele ID (`CA123456`) to a ClinVar Allele ID requires querying the ClinGen Allele Registry REST API (`/allele/CA123456`). This HTTP call is made once per unique ClinGen ID per run, with results cached in:

1. **In-process dict** — lives for the duration of the script run.
2. **Redis** (via `src/lib/clingen.py`) — persistent across runs with the same TTL as other ClinGen allele responses (default 1 day). If Redis is unavailable the in-process cache still prevents duplicate requests within a run.

---

## CLI options

| Option | Default | Description |
|---|---|---|
| `--clinvar-version YYYYMM` | `202601` | ClinVar monthly release to annotate against |
| `--clinvar-namespace NS` | `clinvar` | Prefix for output column names |
| `--cache-dir DIR` | `$CLINVAR_CACHE_DIR` or `/tmp/clinvar_cache` | Directory for downloaded ClinVar TSV files |
| `--dna-clingen-allele-id-col COL` | `dna_clingen_allele_id` | Input column containing ClinGen allele IDs |
| `--delimiter CHAR` | `\t` (TAB) | Field delimiter for input and output |
| `--max-workers N` | `8` | Concurrent worker threads |
| `--skip N` | `0` | Skip the first N data rows |
| `--limit N` | no limit | Stop after processing N rows |
| `--log-level` | `INFO` | Logging verbosity (`DEBUG`, `INFO`, `WARNING`, `ERROR`) |
| `--csv-field-size-limit BYTES` | system default | Increase for large HGVS fields |

---

## Usage

```bash
src/scripts/run_annotate_clinvar.sh input.tsv output.tsv \
    --clinvar-version 202601 \
    --cache-dir ./clinvar_cache
```

To annotate against multiple ClinVar releases, run the script twice with different `--clinvar-version` and `--clinvar-namespace` values:

```bash
src/scripts/run_annotate_clinvar.sh input.tsv output_v1.tsv \
    --clinvar-version 202601

src/scripts/run_annotate_clinvar.sh output_v1.tsv output_v2.tsv \
    --clinvar-version 202504
```

---

## Troubleshooting

**All annotation columns are empty despite valid `dna_clingen_allele_id`**

- Check that the ClinGen API is reachable (`reg.genome.network`).
- Some ClinGen alleles have no corresponding ClinVar record; this is expected for novel variants.
- If many rows are empty, check that `dna_clingen_allele_id` is CA-prefixed and not blank.

**Download fails for a specific `--clinvar-version`**

- NCBI's archive may not yet have a file for that month. Try the most recent published month.
- Confirm internet access from the Docker container.

**`--clinvar-version` format error**

The value must be exactly 6 digits in YYYYMM format (e.g. `202601`). The month must be 01–12.

**Cache directory permission error**

Ensure the `--cache-dir` path is writable from the Docker container. In Docker Compose the cache is mapped to a named volume; run `docker volume ls | grep clinvar` to confirm it exists.
