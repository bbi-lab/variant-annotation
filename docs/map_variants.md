# map_variants

Maps input variants to human-genome reference HGVS strings on GRCh38 and populates ClinGen Allele Registry metadata.

This is the first step in the annotation pipeline. Every downstream step depends on the columns it produces.

---

## Overview

`map_variants` processes each input row according to the type of HGVS data it contains (see [Variant cases](#variant-cases)) and, for each row, does two things:

1. **Derives an assay-level genomic HGVS string** — either by normalising a supplied transcript HGVS through the ClinGen Allele Registry or by aligning a target sequence to GRCh38 via BLAT + VRS mapping (powered by `dcd_mapping`).
2. **Queries the ClinGen Allele Registry** with the assay-level HGVS to retrieve canonical `mapped_hgvs_g`, `mapped_hgvs_c`, and `mapped_hgvs_p` strings, plus the ClinGen allele ID.

---

## Usage

```bash
# Via Docker (recommended):
src/scripts/run_map_variants.sh <input-file> <output-file> [options...]

# Directly inside the container:
python -m src.map_variants <input-file> <output-file> [options...]
```

File paths passed to `run_map_variants.sh` are automatically translated to container-internal paths (`/work/` for data files, `/usr/src/app/` for source-tree files).

---

## Variant cases

The script detects three mutually exclusive cases from the `raw_hgvs_nt` and `raw_hgvs_pro` columns:

### Case 1 — Reference-based transcript HGVS

**Condition:** `raw_hgvs_nt` contains a transcript accession prefix (e.g. `NM_000277.3:c.1218G>A` or `ENST00000316054.9:c.1142G>A`).

**Processing:**
1. The HGVS is passed through `dcd_mapping.vrs_map.fetch_clingen_genomic_hgvs` to obtain an assay-level genomic HGVS (this normalises the notation and resolves any Ensembl → RefSeq conversion if needed).
2. The assay-level HGVS is queried in the ClinGen Allele Registry.
3. The transcript allele matching the original accession is selected from the response to populate `mapped_hgvs_c` and `mapped_hgvs_p`; the GRCh38 genomic allele populates `mapped_hgvs_g`.

**Requirements:** None beyond the input column. `target_sequence` is not required.

**Note:** If `dcd_mapping` is unavailable (e.g. not installed), the script falls back to querying ClinGen directly with the raw input HGVS.

---

### Case 2 — Sequence-based nucleotide HGVS (no accession prefix)

**Condition:** `raw_hgvs_nt` is present but contains no colon (e.g. `c.1218G>A`).

**Processing:**
1. All rows sharing the same `--group-by` column value are batched together for a single BLAT alignment run.
2. `dcd_mapping` aligns the `target_sequence` to GRCh38, selects transcripts, and performs VRS mapping for every row in the batch.
3. The assay-level HGVS from VRS mapping is queried in ClinGen to populate all three output columns.
4. Multi-variant intra-codon haplotypes (e.g. `c.[1A>G;3G>T]`) are normalised to a `delins` expression before mapping when all components fall within one codon.

**Requirements:** `target_sequence` column (or `--targets-file`). `dcd_mapping` must be installed with its data dependencies.

---

### Case 3 — Protein-only variant

**Condition:** `raw_hgvs_nt` is absent (or a blank sentinel such as `_wt`, `_sy`, `=`), and `raw_hgvs_pro` is present (e.g. `p.Ala406Thr` or `p.A406T`).

**Processing:**
1. Protein HGVS strings in one-letter amino acid code (`p.A406T`) are automatically converted to three-letter code (`p.Ala406Thr`) before processing.
2. `dcd_mapping` runs at the protein annotation layer to map the variant.
3. ClinGen is queried to populate `mapped_hgvs_p`. `mapped_hgvs_c` and `mapped_hgvs_g` will be empty unless a protein allele entry in ClinGen carries them (uncommon).

**Requirements:** `target_sequence` column (or `--targets-file`). `dcd_mapping` must be installed.

---

## Input

| Column | Required | Description |
|---|---|---|
| `raw_hgvs_nt` | Conditional | Raw nucleotide HGVS (determines case 1 or 2). Blank sentinels: ``, `NA`, `N/A`, `None`, `_wt`, `_sy`, `=` |
| `raw_hgvs_pro` | Conditional | Raw protein HGVS (used when `raw_hgvs_nt` is absent — case 3). |
| `target_sequence` | Cases 2 & 3 | Full nucleotide or amino acid target sequence used for BLAT alignment. Not needed for case 1. May be supplied via `--targets-file` instead. |
| `target_name` | With `--targets-file` | Key column used to join input rows to the targets file. |

At least one of `raw_hgvs_nt` or `raw_hgvs_pro` must be non-blank. Rows where both are blank receive `mapping_error = "No usable HGVS variant data in this row"`.

---

## Output columns

The following columns are appended to every input row:

| Column | Description |
|---|---|
| `mapped_hgvs_g` | GRCh38 genomic HGVS (e.g. `NC_000012.12:g.102917016C>A`) |
| `mapped_hgvs_c` | Transcript HGVS from ClinGen (e.g. `NM_000277.3:c.1218G>A`). May be empty for protein-only rows. |
| `mapped_hgvs_p` | Protein HGVS from ClinGen (e.g. `NP_000268.1:p.Arg406Gln`). May be empty for non-coding or purely genomic variants. |
| `clingen_allele_id` | ClinGen Allele Registry identifier (e.g. `CA123456`). Empty when ClinGen has no record. |
| `mapping_error` | Non-empty if the row could not be mapped. Downstream steps should check this column. |
| `mapping_warnings` | Non-fatal warnings from HGVS extraction (e.g. coordinate parsing issues). |
| `dna_vrs_digest` | VRS digest for the DNA/genomic allele (populated for cases 2 and 3). |
| `protein_vrs_digest` | VRS digest for the protein allele (populated for case 3). |
| `strand` | Chromosomal strand: `1` (plus) or `-1` (minus). Derived from UTA for case 1, from the BLAT alignment for cases 2 and 3. |

Columns specified in `--drop-columns` are excluded from the output entirely. In the pipeline, `target_sequence` is typically dropped to avoid carrying large sequence blobs through downstream steps.

---

## Options

### Input column names

| Option | Default | Description |
|---|---|---|
| `--raw-hgvs-nt COLUMN` | `raw_hgvs_nt` | Column containing raw nucleotide HGVS strings |
| `--raw-hgvs-pro COLUMN` | `raw_hgvs_pro` | Column containing raw protein HGVS strings |
| `--target-sequence COLUMN` | `target_sequence` | Column containing target sequences |
| `--target-type COLUMN` | `target_type` | Column distinguishing reference from target-sequence rows (informational only) |

### Output column names

| Option | Default | Description |
|---|---|---|
| `--mapped-hgvs-g COLUMN` | `mapped_hgvs_g` | Output column for genomic HGVS |
| `--mapped-hgvs-c COLUMN` | `mapped_hgvs_c` | Output column for transcript HGVS |
| `--mapped-hgvs-p COLUMN` | `mapped_hgvs_p` | Output column for protein HGVS |
| `--mapping-error COLUMN` | `mapping_error` | Output column for error messages |
| `--mapping-warnings COLUMN` | `mapping_warnings` | Output column for non-fatal warnings |
| `--dna-vrs-digest COLUMN` | `dna_vrs_digest` | Output column for DNA VRS digests |
| `--protein-vrs-digest COLUMN` | `protein_vrs_digest` | Output column for protein VRS digests |
| `--strand COLUMN` | `strand` | Output column for chromosomal strand |
| `--clingen-allele-id COLUMN` | `clingen_allele_id` | Output column for ClinGen allele IDs |

### Grouping and targets

| Option | Default | Description |
|---|---|---|
| `--group-by COLUMN` | `target_sequence` | Column used to batch sequence-based rows (cases 2 & 3). Rows sharing the same value are aligned together in a single BLAT run. Use `gene_symbol` or a dataset name rather than `target_sequence` when multiple rows share the same sequence, to avoid redundant alignments. |
| `--targets-file FILE` | — | Optional TSV/CSV file containing target sequences and other target-level columns. Joined to the input on the column specified by `--target-name`. Columns from the targets file are merged onto the input row before processing; input columns take precedence when both are present. |
| `--target-name COLUMN` | `target_name` | Join column for `--targets-file`. |

### Row selection

| Option | Default | Description |
|---|---|---|
| `--skip N` | `0` | Skip the first N data rows. Use to resume an interrupted run. |
| `--limit N` | — | Process at most N data rows after `--skip`. |

### Concurrency and retries

| Option | Default | Description |
|---|---|---|
| `--max-clingen-concurrency N` | `5` | Maximum concurrent ClinGen Allele Registry HTTP requests. |
| `--dcd-chunk-on-137 / --no-dcd-chunk-on-137` | on | Automatically retry groups that fail with BLAT exit code 137 (OOM) using progressively smaller chunks. |
| `--dcd-chunk-size-on-137 N` | `500` | Initial chunk size for the first retry on BLAT error 137. |
| `--dcd-max-retry-attempts N` | `3` | Maximum retry attempts before the group is marked as failed. |

### Output ordering

| Option | Default | Description |
|---|---|---|
| `--preserve-order MODE` | `groups` | `groups` — preserve input order while processing groups (recommended). `index` — buffer all results by row index and emit in strict input order. `no` — write each result immediately as it arrives (may emit out of input order). |

### Column management

| Option | Default | Description |
|---|---|---|
| `--drop-columns COLUMN` | — | Omit a column from the output. Repeat to drop multiple columns. In the pipeline, `target_sequence` is routinely dropped. |

### Reusing prior results

| Option | Default | Description |
|---|---|---|
| `--merge-existing FILE` | — | Reuse mapped results from a previously annotated file. Rows in the input are matched by `(raw_hgvs_nt, raw_hgvs_pro)` (plus any `--merge-match-col` columns). Matching rows are written to the output using the prior result and are skipped from fresh processing. Repeat to supply multiple merge sources; the first match wins. |
| `--merge-match-col COLUMN` | — | Additional column that must match for merge reuse. Repeat for multiple columns. |

### Shell wrapper extras

These options are consumed by `run_map_variants.sh` and are not passed to the Python script:

| Option | Description |
|---|---|
| `--rebuild-image` | Rebuild the Docker image before running. |
| `--no-build-cache` | Pass `--no-cache` to `docker compose build` (use with `--rebuild-image`). |

---

## ClinGen query batching

ClinGen Allele Registry queries are issued concurrently within each processing
unit to reduce wall-clock time.

| Case | What is batched |
|---|---|
| **1** (reference-based) | After `fetch_clingen_genomic_hgvs` resolves each assay-level HGVS (sequentially — this is a synchronous `dcd_mapping` call), all unique assay-level strings in the input chunk are collected and queried concurrently. |
| **2** (sequence-based) | After `dcd_mapping` completes for the whole group, all unique assay-level strings are collected and queried concurrently. |
| **3** (protein-only) | Same as case 2 — batched after the group's DCD pipeline run finishes. |

The `--max-clingen-concurrency` option (default: 5) sets the `asyncio.Semaphore`
limit — the maximum number of in-flight HTTP requests at any moment. Lower values
reduce the chance of hitting ClinGen rate limits; higher values reduce wall time
on fast networks.

Debug logs show per-group batch timing when `--log-level DEBUG` is set:

```
DEBUG: Batch querying ClinGen for 47 variants in group 'TP53' (transcript: NM_000546.6)
DEBUG: ClinGen batch query completed for 47 variants in 9.45 seconds
```

---

## Dependencies

### Always required
- **ClinGen Allele Registry** — public REST API at `https://reg.clinicalgenome.org`. Internet access required. Requests are retried automatically with backoff.
- **Redis** (optional) — if a Redis service is reachable, ClinGen responses and intermediate genomic HGVS lookups are cached for the duration of `VEP_CACHE_TTL_SECONDS` (or the ClinGen equivalent), avoiding redundant API calls across runs.

### Required for cases 2 and 3 only
- **`dcd_mapping`** — must be installed in the container image. Install from the `mavedb-main` branch:
  ```
  pip install 'dcd-mapping @ git+https://github.com/VariantEffect/dcd_mapping2.git@mavedb-main'
  ```
- **SeqRepo snapshot** — mounted at `/usr/local/share/seqrepo` inside the container. Managed via the `seqrepo` Docker service.
- **UTA database** — PostgreSQL database accessible via `UTA_DB_URL`. Managed via the `uta` Docker service.
- **Gene Normalizer database** — used by `cool_seq_tool` (a `dcd_mapping` dependency). Managed via the `db` Docker service.
- **BLAT binary** — UCSC BLAT tool. Included in the `map-variants` Docker image. The GRCh38 2-bit genome file is downloaded automatically to the `dcd_mapping` local store on first use.

> **Apple Silicon note:** The `map-variants` service runs as `linux/amd64` (via Rosetta emulation) because BLAT does not have an ARM64 build. This is handled automatically by the Docker Compose profile.

---

## Environment variables

| Variable | Description |
|---|---|
| `UTA_DB_URL` | PostgreSQL connection URL for the UTA transcript database (e.g. `postgresql://uta_admin@uta/uta/uta_20241220`). |
| `SEQREPO_DATA_PROXY_URL` | SeqRepo access URL for `cool_seq_tool`. |
| `GENE_NORM_DB_URL` | PostgreSQL URL for the Gene Normalizer database. |
| `REDIS_URL` | Redis connection URL. Used by the ClinGen fetch cache patch to avoid redundant HTTP calls. |

These are typically set via `settings/.env` and injected automatically by Docker Compose.

---

## Typical pipeline invocation

From `data/cvfg/generate.sh` (Step 1):

```bash
src/scripts/run_map_variants.sh data/cvfg/v7/cvfg_variants.0.tsv data/cvfg/v7/cvfg_variants.1b.tsv \
  --drop-columns target_sequence \
  --max-clingen-concurrency 3 \
  --skip 90947
```

- `--drop-columns target_sequence` removes the large sequence column from the output to keep file sizes manageable.
- `--max-clingen-concurrency 3` limits simultaneous ClinGen requests (the default of 5 can trigger rate-limit errors for large inputs).
- `--skip 90947` resumes from a row where a previous run was interrupted.

When running for the first time against a new input file:

```bash
src/scripts/run_map_variants.sh variants.0.tsv variants.1.tsv \
  --drop-columns target_sequence
```

---

## Resuming an interrupted run

The script processes rows in a single pass and writes output incrementally (with `--preserve-order groups`). If a run is interrupted:

1. Check how many rows the output file contains (`wc -l output.tsv` minus 1 for the header).
2. Restart with `--skip <rows_already_written>` and append to a new file, then concatenate the header + the new output rows.

Alternatively, use `--merge-existing` on a subsequent full run: matched rows are taken from the prior output, and only unprocessed rows are sent to `dcd_mapping`/ClinGen.

```bash
src/scripts/run_map_variants.sh variants.0.tsv variants.1.tsv \
  --drop-columns target_sequence \
  --merge-existing variants.1.partial.tsv
```

---

## BLAT exit code 137 (OOM)

BLAT is invoked as a subprocess and can be killed by the kernel OOM killer when the reference genome is too large to fit in memory. Symptoms:

- `BLAT alignment failed` errors in the log
- Error message containing `error code 137`

Mitigations (applied automatically):

1. `--dcd-chunk-on-137` (on by default) splits the failing group into smaller batches and retries up to `--dcd-max-retry-attempts` times with halving chunk sizes.
2. Increase Docker Desktop memory allocation (Settings → Resources). 8 GB minimum; 16 GB recommended for large gene groups.

If chunking still fails, the affected group is recorded in `mapping_error` and processing continues.

---

## Troubleshooting

**`dcd_mapping` not found**
: Cases 2 and 3 require `dcd_mapping`. The map-variants Docker image includes it. If running outside Docker, install it manually (see [Dependencies](#dependencies)).

**`target_sequence` missing for case 2/3 row**
: The row will receive `mapping_error = "target_sequence is required for case N but column 'target_sequence' is missing or blank"`. Supply the sequence directly in the input file or via `--targets-file`.

**ClinGen returned no data**
: Transient network errors are retried with exponential backoff. If ClinGen consistently returns nothing for a specific HGVS, the variant may not be registered; check the ClinGen Allele Registry directly.

**Ensembl ENST accessions**
: If `raw_hgvs_nt` contains an ENST-prefixed HGVS (e.g. `ENST00000316054.9:c.1142G>A`), the script attempts to resolve a RefSeq NM_ accession via the Ensembl REST API before querying ClinGen. This requires internet access to `https://rest.ensembl.org`.

**All mapping columns blank but no error recorded**
: This is caught as `mapping_error = "Mapping produced no HGVS output"`. It can happen when a ClinGen CA allele record exists but contains no GRCh38 genomic coordinates and no transcript alleles matching the input accession.
