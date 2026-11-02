# fetch_mavedb_scores

`src/fetch_mavedb_scores.py` builds the pipeline's starting variant table from a list of MaveDB score sets, using the public MaveDB REST API. The output goes straight into `map_variants`. Its columns are those of `cvfg_variants.0.tsv` (exported from a MaveDB database copy by `data/cvfg/mavedb_queries.sql`), minus the ones the pipeline doesn't use; see [Differences from the database export](#differences-from-the-database-export).

No dataset-specific handling is applied. If a score set needs a corrected target sequence, deduplication, or edited HGVS strings, apply that to the output before mapping.

## Usage

```bash
src/scripts/run_fetch_mavedb_scores.sh data/excalibr_domainome/excalibr_domainome.0.tsv \
    --score-sets-file data/excalibr_domainome/score_sets.tsv

# Or name score sets directly (repeatable; combined with --score-sets-file if both are given):
src/scripts/run_fetch_mavedb_scores.sh data/test/variants.0.tsv \
    --score-set-urn urn:mavedb:00000895-a-1 \
    --score-set-urn urn:mavedb:00000736-a-1
```

Without Docker: `python -m src.fetch_mavedb_scores <output-file> [options]`.

## Score sets file

A CSV/TSV (delimiter from the extension) with a `score_set_urn` column. Other columns are ignored, so the same file can carry `dataset_name` and `requested_calibration_urn` for later steps (`annotate_mavedb --requested-calibrations-file`).

If the file has a `preferred_transcript` column (name set by `--preferred-transcript-col`), its value is copied to every variant of that score set. Pass it on to `map_variants --preferred-transcript-col preferred_transcript` to fix the transcript used for mapped `c.`/`p.` strings, instead of relying on automatic MANE selection.

```
dataset_name	score_set_urn	requested_calibration_urn	preferred_transcript
MECP2_Beltran_2024_MBD	urn:mavedb:00000895-a-1		NM_001110792.2
```

## What is fetched

For each score set:

| Request | Used for |
|---|---|
| `GET /api/v1/score-sets/{urn}` | `score_set_title`, `reference_based`, `target_sequence_type`, `target_sequence` (first target gene only; a warning is logged if there are more) |
| `GET /api/v1/score-sets/{urn}/scores` (CSV) | `variant_urn`, `raw_hgvs_nt`, `raw_hgvs_pro`, `score`, `rna_score*`, `score_data` |
| `GET /api/v1/score-sets/{urn}/variant-details` (NDJSON), only with `--include-mavedb-mappings` | `mavedb_mapped_hgvs_g`/`_c`/`_p` from MaveDB's own mapping |

plus one `GET /api/v1/score-sets/mapped-genes` call for `gene_symbol`. If MaveDB has no mapped gene for a score set, the target gene's `mappedHgncName` is used, and failing that `gene_symbol` is left blank with a warning. The target gene's display name (e.g. "MECP2 Methyl-CpG DNA binding domain") is never used as a gene symbol.

- `score` and `rna_score*` keep MaveDB's text exactly; `NA` becomes blank.
- `score_data` is a JSON object of every `/scores` column other than `accession`, `hgvs_nt`, `hgvs_splice`, and `hgvs_pro`, with numbers as numbers and `NA` values omitted. `rna_score` is taken from a `rna_score` or `score_rna` column.
- Rows are ordered by gene symbol, then by score set (in input order), then by variant number.

## Output columns

`gene_symbol`, `score_set_title`, `variant_urn`, `raw_hgvs_nt`, `raw_hgvs_pro`, `score`, `rna_score`, `rna_score_d6`, `rna_score_d20`, `score_data`, `reference_based`, `target_sequence_type`, `target_sequence`, `preferred_transcript`.

With `--include-mavedb-mappings`, `mavedb_mapped_hgvs_g`, `mavedb_mapped_hgvs_c`, and `mavedb_mapped_hgvs_p` follow `reference_based`. These hold MaveDB's own mapping of each variant, for reference; the pipeline doesn't use them, since `map_variants` computes its own `mapped_hgvs_*` from `raw_hgvs_nt`/`raw_hgvs_pro` and `target_sequence`. Variants MaveDB could not map (typically `p.=`/`_wt` rows) are absent from `/variant-details` and get blank values, and a score set with no mapped variants (HTTP 404) gets all three blank. Fetching `/variant-details` is the largest request per score set, so it is skipped unless the option is given.

## Differences from the database export

- `variant_id` (MaveDB's internal integer ID) and `mavedb_mapping_error` are omitted. The API exposes neither, and the pipeline doesn't use them.
- `mavedb_mapped_hgvs_g`/`_c`/`_p` are only included with `--include-mavedb-mappings`.

## Options

| Option | Default | Description |
|---|---|---|
| `--score-sets-file PATH` | — | CSV/TSV with a `score_set_urn` column |
| `--score-set-urn URN` | — | Score set to fetch; repeatable |
| `--preferred-transcript-col COLUMN` | `preferred_transcript` | Score sets file column copied to the output's `preferred_transcript` |
| `--include-mavedb-mappings` | off | Add `mavedb_mapped_hgvs_g`/`_c`/`_p` from MaveDB's own mapping |
| `--mavedb-api-url URL` | `MAVEDB_API_URL` env var, or `https://api.mavedb.org` | MaveDB REST API base URL |
| `--retries N` | `8` | Retries per request on HTTP 500/502/503/504, with exponential backoff |
| `--timeout SECONDS` | `300` | Per-request timeout |
| `--sleep-between-requests SECONDS` | `0` | Pause between score sets |
| `--log-level LEVEL` | `INFO` | Logging verbosity |

At least one of `--score-sets-file` or `--score-set-urn` is required.
