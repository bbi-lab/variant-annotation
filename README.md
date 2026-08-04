# variant-annotation

## License

This project is licensed under the GNU Affero General Public License v3.0 or
later. See the [LICENSE](LICENSE) file.

The AGPL choice reflects that parts of the variant-mapping pipeline are
intentionally modeled on established MaveDB workflow behavior and may include
adapted implementation ideas or logic in that area of the codebase.

## Prerequisites

- [Docker](https://docs.docker.com/get-docker/) and [Docker Compose](https://docs.docker.com/compose/install/)
- Python 3.11+ (for local development; Docker provides 3.11 for pipeline scripts)
- Git
- ~20 GB free disk space (for databases, caches, and resources)

## Development setup

All pipeline scripts and services run in Docker. The development environment is defined in `compose.yaml` and includes the following services:

**Annotation pipeline tools** (run with `--profile tools`):
- **map-variants**: Maps variants to GRCh38 HGVS (includes BLAT, SeqRepo, Gene Normalizer)
- **reverse-translate-proteins**: Reverse-translates protein variants to DNA candidates
- **annotate-clinvar**: Annotates with ClinVar clinical significance
- **annotate-gnomad**: Annotates with gnomAD allele frequencies
- **annotate-spliceai**: Annotates with SpliceAI splice-impact scores
- **annotate-vep**: Annotates with VEP mutational consequence
- **flatten-dna-variants**: Flattens multi-candidate DNA variants to one row per variant

**Data and dependency services** (always running):
- **cdot**: Transcript coordinate data service
- **db**: PostgreSQL database (Gene Normalizer backend)
- **redis**: Cache for cdot service
- **seqrepo**: Sequence repository (genome references)
- **uta**: UTA database (transcript/protein lookups)

### 0. Clone the repository

```bash
git clone https://github.com/bbi-lab/variant-annotation.git
cd variant-annotation
```

### 1. Configure environment variables

Copy the template environment file:

```bash
cp settings/.env.template settings/.env
```

The template includes sensible defaults for development (including `VARIANT_DATA_DIR=./data`). You can optionally customize values in `.env` if needed:

```bash
# Override the data directory for TSV input/output (default: ./data)
VARIANT_DATA_DIR=./data

# Optional: Set custom UTA database URL
# UTA_DB_URL=postgresql://user:password@host:port/uta

# Optional: Postgres credentials (defaults are fine for local dev)
# POSTGRES_USER=postgres
# POSTGRES_PASSWORD=postgres
```

### 2. Prepare the UTA database (one-time)

Before first startup, download the UTA database dump (~500 MB):

```bash
src/scripts/fetch_uta_dump.sh
```

This step is required for transcript lookups in `reverse_translate_protein_variants` and `add_vcf_identifiers`. It downloads the dump to a volume-mounted location in the Docker container.

### 3. Build and start services

```bash
docker compose up --build
```

This will:
- Build the Docker image (including BLAT, samtools, htslib, Java for Hail)
- Start all services (containers for app, map-variants, uta, seqrepo, postgres)
- Initialize databases on first run

**Expected output**: Services will be ready when you see messages indicating database initialization is complete. This may take 5-10 minutes on first run.

The `cdot` service is built from source at [github.com/VariantEffect/cdot_rest](https://github.com/VariantEffect/cdot_rest) — no manual clone required.

### 4. Open a shell in the app container

In a new terminal (while services are running):

```bash
docker compose exec app bash
```

This opens an interactive shell in the app container where you can run Python commands, tests, or scripts. Your code changes in the host repository are immediately visible in the container (via bind mount).

### 5. Run pipeline scripts on a TSV file

Put your input file in `./data` (or your `VARIANT_DATA_DIR` folder), then run:

```bash
src/scripts/run_map_variants.sh input.tsv output.tsv
```

**File mounting:**
- If `input.tsv` exists in the repository working tree, the wrapper script automatically uses the project bind mount (`/usr/src/app/...`)
- Otherwise it defaults to the `/work` mount backed by `./data` (or `VARIANT_DATA_DIR`)
- Output files are written back to the same host-mounted folder

**Platform-specific notes:**
- On Apple Silicon, `map-variants` runs as `linux/amd64` (via emulation) so UCSC BLAT works
- The Docker image installs `polars[rtcompat]` to avoid AVX-related crashes under emulation
- Class-1 (accession-referenced) nucleotide HGVS variants are normalized through `dcd_mapping` (matching MaveDB behavior) before ClinGen extraction, which preserves accession/version handling for ENST inputs

**Example with options:**
```bash
src/scripts/run_map_variants.sh input.tsv output.tsv \
    --group-by gene_symbol \
    --raw-hgvs-nt raw_hgvs_nt \
    --raw-hgvs-pro raw_hgvs_pro
```

**Reuse existing mapped results** to avoid reprocessing matching variants:

```bash
src/scripts/run_map_variants.sh input.tsv output.tsv \
	--merge-existing prior_results_a.tsv \
	--merge-existing prior_results_b.tsv \
	--merge-match-col gene_symbol
```

Rows are merged only when `raw_hgvs_nt`, `raw_hgvs_pro`, and any
`--merge-match-col` columns are identical between input and merge source rows.
Other columns (including `target_sequence`) are ignored for merge matching.
If a match is found, the previously mapped values are copied directly to output and
the row is skipped from new mapping work.

The output file is written back to the same host-mounted folder.

### BLAT Error 137 Retry Strategy

When processing large groups of variants, the BLAT process may fail with error code 137
(SIGKILL), typically caused by insufficient memory. The `map-variants` command includes
automatic error-driven retry with progressive chunking to handle this gracefully.

Retry options:

- `--dcd-chunk-on-137` / `--no-dcd-chunk-on-137` (default: enabled)
  Automatically retry groups that fail with BLAT error 137 by splitting them into
  progressively smaller chunks.

- `--dcd-chunk-size-on-137 SIZE` (default: 500)
  Initial chunk size for the first retry attempt. On subsequent retries, the chunk
  size is halved each time (500 → 250 → 125, etc.).

- `--dcd-max-retry-attempts N` (default: 3)
  Maximum number of retry attempts allowed before recording the group as failed.

Example with retry tuning:

```bash
src/scripts/run_map_variants.sh input.tsv output.tsv \
	--dcd-chunk-on-137 \
	--dcd-chunk-size-on-137 300 \
	--dcd-max-retry-attempts 4
```

For more details on error handling, see the **Troubleshooting** section below.

## Annotation Pipeline

The annotation pipeline comprises multiple stages that build upon each other. All scripts process Tab-Separated Values (TSV) or Comma-Separated Values (CSV) files. Each stage will be run using wrapper scripts in `src/scripts/`.

### Pipeline Overview

The complete pipeline flow is:

```
Input variants
    ↓
[1] map_variants ──→ mapped HGVS (c./g./p.)
    ↓
[2] reverse_translate_protein_variants ──→ protein candidates in c./g. format
    ↓
[3] add_dna_clingen_allele_ids ──→ DNA-level ClinGen allele IDs
    ↓
[4] add_vcf_identifiers (optional) ──→ parsed positions/alleles
    ↓
[5] annotate_clinvar (optional) ──→ ClinVar clinical significance
    ↓
[6] annotate_gnomad (optional) ──→ gnomAD allele frequencies
    ↓
[7] annotate_spliceai (optional) ──→ SpliceAI delta scores
    ↓
[8] annotate_erepo (optional) ──→ expert-panel classifications
    ↓
[9] annotate_vep (optional) ──→ VEP mutational consequence
    ↓
[10] annotate_mavedb (optional) ──→ MaveDB functional classifications
    ↓
[11] annotate_predictors (optional) ──→ REVEL / AlphaMissense / MutPred2 scores
    ↓
[12] flatten_dna_variants (optional) ──→ flattened DNA-only variants
    ↓
Output: one row per DNA variant (or fully annotated variants)
```

For a detailed comparison with MaveDB's annotation workflow, see [docs/differences_from_mavedb.md](docs/differences_from_mavedb.md).

### Step 1: Map Variants (Required)

**Purpose:** Normalize input variants (nucleotide or protein HGVS) to human-genome reference HGVS strings on GRCh38 and populate ClinGen Allele Registry metadata.

See [docs/map_variants.md](docs/map_variants.md) for full reference documentation including all three input cases, CLI options, and troubleshooting.

**Input columns:** `raw_hgvs_nt` (nucleotide HGVS), `raw_hgvs_pro` (protein HGVS), `target_sequence` (required for sequence-based and protein-only rows)

**Output columns:** `mapped_hgvs_g`, `mapped_hgvs_c`, `mapped_hgvs_p`, `clingen_allele_id`, `mapping_error`, `mapping_warnings`, `dna_vrs_digest`, `protein_vrs_digest`, `strand`

**Command:**
```bash
src/scripts/run_map_variants.sh input.tsv output.tsv \
    --drop-columns target_sequence
```

**Notes:**
- Three variant cases are detected automatically from the input columns — see [docs/map_variants.md](docs/map_variants.md) for full details:
  - **Case 1** (`raw_hgvs_nt` with a transcript accession prefix, e.g. `NM_000277.3:c.1218G>A`): normalised through `dcd_mapping` and ClinGen; no target sequence needed
  - **Case 2** (`raw_hgvs_nt` without an accession prefix, e.g. `c.1218G>A`): BLAT-aligned from `target_sequence`; all rows sharing the same `--group-by` value are aligned together
  - **Case 3** (protein-only, `raw_hgvs_pro` with no `raw_hgvs_nt`): BLAT-aligned at the protein annotation layer
- `--drop-columns target_sequence` removes the large sequence column from the output (recommended)
- Use `--group-by gene_symbol` (or another stable column) instead of the default `target_sequence` when many rows share the same sequence — this avoids redundant BLAT alignments
- Use `--max-clingen-concurrency 3` to reduce the chance of rate-limit errors from the ClinGen API for large inputs (default: 5)
- Use `--skip N` to resume an interrupted run from row N
- Use `--merge-existing prior_output.tsv` to reuse results from a previous partial run without re-processing matched rows
- Use `--preferred-transcript NM_ACCESSION` (e.g. `--preferred-transcript NM_007194.4`) to force a specific MANE Select transcript for all groups when automatic selection picks the wrong one
- Use `--preferred-transcript-col COLUMN` to specify the preferred transcript per row from a column in the input file; blank values fall back to automatic selection (or `--preferred-transcript` if also provided)
- See [docs/map_variants.md](docs/map_variants.md) for all options, dependency setup, and troubleshooting
- See [BLAT Error 137 Retry Strategy](#blat-error-137-retry-strategy) for handling memory issues

### Step 2: Reverse-Translate Protein Variants (Conditional)

**Purpose:** For protein-only variants mapped in step 1, reverse-translate them into every DNA (c./g.) codon substitution that produces the observed amino acid change. A single protein change can arise from two or three synonymous codons; this step enumerates all of them so downstream annotation (ClinVar, gnomAD, VEP) has a DNA variant to query.

See [docs/reverse_translate_protein_variants.md](docs/reverse_translate_protein_variants.md) for full reference documentation including transcript resolution order, batching strategy, and troubleshooting.

**Input columns:** `mapped_hgvs_p` (non-empty), `mapped_hgvs_c` and `mapped_hgvs_g` (must be empty/blank for this step to apply)

**Output columns:** Updates `mapped_hgvs_c` and `mapped_hgvs_g` with pipe-delimited candidates (e.g., `NM_000277.3:c.1216G>A|NM_000277.3:c.1217C>A`). Also adds `assayed_variant_level`, `reverse_translation_error`, `reverse_translation_warnings`, and parsed position/allele columns for each candidate (`mapped_hgvs_c_start`, `mapped_hgvs_g_ref`, etc.)

**Command:**
```bash
src/scripts/run_reverse_translate_protein_variants.sh output.tsv output_rt.tsv
```

**Notes:**
- Only processes rows where `mapped_hgvs_p` is non-empty and both `mapped_hgvs_c` and `mapped_hgvs_g` are blank
- Pipe-delimited candidates are position-aligned across all output columns: an empty slot (e.g. `NM_...:c.1216G>A||NM_...:c.1218G>A`) means the middle candidate could not be resolved
- DNA variants from step 1 pass through unchanged; their parsed position columns are also populated
- Requires the `reverse-translate-variants` CLI (installed in the Docker image) and UTA database access
- Use `--transcript-fallback-column raw_hgvs_nt` if the protein HGVS string lacks a transcript accession prefix
- Use `--include-indels` to also generate small insertion/deletion candidates (off by default)
- Use `--wt-codon-mode unambiguous` or `--wt-codon-mode all` (requires `--include-indels`) to append WT codon candidates for synonymous protein variants (`ref AA == alt AA`, e.g. `p.Met1Met`)
- See [docs/reverse_translate_protein_variants.md](docs/reverse_translate_protein_variants.md) for all options, transcript resolution logic, and troubleshooting

### Step 3: Add DNA-Level ClinGen Allele IDs (Required for annotation)

**Purpose:** Resolve a DNA-level ClinGen Allele Registry (CAR) identifier for every variant. For protein-origin rows with multiple reverse-translated candidates, each candidate gets its own ID. The resulting `dna_clingen_allele_id` column is the primary key for all downstream annotation steps.

See [docs/add_dna_clingen_allele_ids.md](docs/add_dna_clingen_allele_ids.md) for full reference documentation including the lookup strategy, caching details, and known-misses file generation.

**Input columns:** `mapped_hgvs_c`, `mapped_hgvs_g`, `clingen_allele_id` (from step 1), `raw_hgvs_nt`, `raw_hgvs_pro`

**Output columns:** `dna_clingen_allele_id` (pipe-delimited, aligned to `mapped_hgvs_c`/`_g` candidates)

**Command:**
```bash
src/scripts/run_add_dna_clingen_allele_ids.sh output_rt.tsv output_clingen.tsv
```

**Notes:**
- **DNA rows** (`raw_hgvs_nt` non-blank): if `clingen_allele_id` already holds a `CA`-prefixed ID and the row has a single candidate, it is reused without a network request
- **Protein-origin rows**: ClinGen is queried for each reverse-translated candidate independently using the c. HGVS first, falling back to g.
- Empty slots preserve cardinality: `CA101||CA103` means 3 candidates, the middle one had no ClinGen match
- ClinGen misses (404 or no allele ID) are cached in Redis as sentinels so repeated runs skip them without another HTTP request
- Use `--known-misses-file` to pre-load a list of HGVS strings known to have no ClinGen record, bypassing even the Redis lookup
- See [docs/add_dna_clingen_allele_ids.md](docs/add_dna_clingen_allele_ids.md) for the full lookup strategy, caching details, and how to generate a known-misses file

### Step 4: Add Parsed Position/Allele Columns (Optional)

**Purpose:** Decompose the three mapped HGVS columns into VCF-style fields — chromosome/transcript accession, start, stop, reference allele, and alternate allele — for easier filtering, joining to other datasets, and variant comparison. Deletions and insertions follow the VCF left-anchor convention.

See [docs/add_vcf_identifiers.md](docs/add_vcf_identifiers.md) for full reference documentation including HGVS parsing rules, VCF anchor convention, and all column definitions.

**Input columns:** `mapped_hgvs_g`, `mapped_hgvs_c`, `mapped_hgvs_p`

**Output columns:** For each of the three HGVS columns: `<col>_chromosome` (or `<col>_transcript` for c.), `<col>_start`, `<col>_stop`, `<col>_ref`, `<col>_alt`. For g. and c. columns these are pipe-delimited when the input is pipe-delimited (matching the candidate cardinality from step 2). Also adds `touches_intronic_region` and `spans_intron` as pipe-delimited per-candidate flags (`"true"`/`"false"`) aligned to `mapped_hgvs_c` (empty slot preserved for empty candidates).

**Command:**
```bash
src/scripts/run_add_vcf_identifiers.sh output_clingen.tsv output_parsed.tsv
```

**Notes:**
- If step 2 (`reverse_translate_protein_variants`) was run, derived columns are already present; this step re-computes them (idempotent for DNA rows, with per-candidate intronic flags preserved and aligned)
- If step 2 was skipped (no protein variants), this is the only step that adds derived position columns
- Inversions are expanded to the reverse-complement sequence; protein HGVS is converted to one-letter amino acid codes
- UTA is required for resolving missing ref alleles in `del`/`dup`/`inv` variants; those fields are left blank if UTA is unavailable
- See [docs/add_vcf_identifiers.md](docs/add_vcf_identifiers.md) for full column list, parsing rules, and VCF anchor details

### Step 5: Annotate with ClinVar Data (Optional)

**Purpose:** Look up clinical significance, review status, and star ratings from ClinVar. Resolves each `dna_clingen_allele_id` to a ClinVar Allele ID via the ClinGen API, then looks up that ID in the cached monthly ClinVar variant-summary TSV from NCBI.

See [docs/annotate_clinvar.md](docs/annotate_clinvar.md) for full reference documentation including star-rating table, ClinGen→ClinVar ID resolution, and caching details.

**Input columns:** `dna_clingen_allele_id` (from step 3)

**Output columns:** `clinvar.<YYYYMM>.variation_id`, `.allele_id`, `.clinical_significance`, `.review_status`, `.stars`, `.last_review_date`

**Command:**
```bash
src/scripts/run_annotate_clinvar.sh output_parsed.tsv output_clinvar.tsv \
    --clinvar-version 202601 \
    --cache-dir ./clinvar_cache
```

**Notes:**
- Pipe-delimited `dna_clingen_allele_id` candidates are each annotated independently; output columns are also pipe-delimited
- Downloads and caches the ClinVar TSV on first run; the cached file is reused permanently (monthly archives are immutable)
- `stars` (0–4) is derived from `review_status` following ClinVar's own star-rating definitions
- Default namespace is `clinvar`; use `--clinvar-namespace` to customise (e.g. when running two versions side by side)
- Use `--dna-clingen-allele-id-col` if your input uses a different column name
- See [docs/annotate_clinvar.md](docs/annotate_clinvar.md) for star-rating table, ClinGen→ClinVar resolution, caching details, and all options

### Step 6: Annotate with gnomAD Allele Frequencies (Optional)

**Purpose:** Look up allele frequencies and population-specific metrics from gnomAD using a local Hail table cache.

See [docs/annotate_gnomad.md](docs/annotate_gnomad.md) for full reference documentation including lookup mode details, CAID resolution strategies, and Athena tuning options.

**Input columns:** `dna_clingen_allele_id` (from step 3)

**Output columns (standard):** `gnomad.<VERSION>.minor_allele_frequency`, `.allele_frequency`, `.allele_count`, `.allele_number`, `.faf95_max`, `.faf95_max_ancestry`, `.filters`, `.exome_filters`, `.genome_filters`, `.gene_symbols`

**Optional histogram columns (Hail mode):** When `--age-histograms` or `--allele-balance-histograms` is passed, per-bin frequency columns are added with names like `gnomad.<V>.age_hist_exome_het.bin_1_17.5_25`, `gnomad.<V>.ab_hist_genome_adj.bin_1_0_0.05`, etc. See [docs/annotate_gnomad.md](docs/annotate_gnomad.md) for full column naming details.

**Command (standard annotation):**
```bash
src/scripts/run_annotate_gnomad.sh output_clinvar.tsv output_final.tsv \
    --gnomad-version v4.1 \
    --cache-dir ./gnomad_cache
```

**First-time setup — build cache with all histogram columns included:**
```bash
src/scripts/run_annotate_gnomad.sh /dev/null /dev/null \
    --gnomad-version v4.1 \
    --cache-dir ./gnomad_cache \
    --download-only \
    --refresh-cache \
    --gnomad-ht-uri gs://gcp-public-data--gnomad/release/4.1/ht/joint/gnomad.joint.v4.1.sites.ht \
    --age-histograms exome,genome,joint \
    --allele-balance-histograms exome,genome,joint
```

**Annotate with all histograms enabled:**
```bash
src/scripts/run_annotate_gnomad.sh output_clinvar.tsv output_final.tsv \
    --gnomad-version v4.1 \
    --cache-dir ./gnomad_cache \
    --age-histograms exome,genome,joint \
    --allele-balance-histograms exome,genome,joint
```

**QC filtering at annotation time:**

Two flags control whether variants that failed gnomAD QC are included in the output.
Both flags treat a missing gnomAD match as "no annotation" (columns left empty) regardless of setting.

| Flag | Semantics |
|------|-----------|
| `--require-pass` | Exclude variants whose gnomAD **combined** `filters` field is non-empty (i.e. the variant failed at least one QC filter across all callsets). Equivalent to `FILTER == PASS` in a joint VCF. |
| `--callset-pass-filter none` _(default)_ | No callset-level filtering; all matched records are annotated. |
| `--callset-pass-filter any` | Only annotate variants that passed QC in **at least one** callset (`exome_filters` or `genome_filters` is empty). |
| `--callset-pass-filter all` | Only annotate variants that passed QC in **both** callsets (`exome_filters` and `genome_filters` are both empty). |

`--require-pass` and `--callset-pass-filter` can be combined; a variant must satisfy both criteria to receive annotation.

`--callset-pass-filter any/all` requires per-callset filter data in the gnomAD cache (i.e. the source Hail table must have separate `exome.filters` / `genome.filters` fields, as in the gnomAD v4.1 joint sites table). An error is raised if the cache was built from a table without those fields; in that case use `--require-pass` instead.

**Notes (Hail mode):**
- Requires Java runtime and Hail library (installed via `gnomad` extra)
- First run downloads the source Hail table from GCP public data; this creates a local indexed cache
- Subsequent runs reuse the local cache (much faster)
- Default gnomAD version is `v4.1`; customize with `--gnomad-version` flag
- Custom table URI: `--gnomad-ht-uri gs://your-bucket/path/to/table.ht`
- Output columns are pipe-delimited and candidate-aligned across all annotation fields
- Supports custom DNA ID column via `--dna-clingen-allele-id-col` if needed
- Cache refresh: use `--refresh-cache` flag to re-download the source table
- Histogram columns (`--age-histograms`, `--allele-balance-histograms`) must be requested both when building the cache and when annotating; passing them to an older cache logs a warning and omits those columns
- Lookup results are also cached in Redis between runs when Redis is available (`GNOMAD_CACHE_REDIS_ENABLED`); useful when annotating overlapping variant sets across multiple runs

#### Athena execution mode (alternative to Hail)

Use `--execution-mode athena` to query gnomAD via AWS Athena instead of a local Hail cache. This is intended for use in environments where the gnomAD data is already registered as an Athena table — for example, the same parquet-backed gnomAD dataset used by a **MaveDB worker** environment. (See the [MaveDB API documentation](https://github.com/VariantEffect/mavedb-api/wiki/gnomAD-Data-Ingestion).) No local Hail installation or GCP access is required.

The Athena table schema and query structure match the `gnomad_variant_data_for_caids` function in `mavedb-api/src/mavedb/lib/gnomad.py`, which selects `caid`, `joint.freq.all.ac`, `joint.freq.all.an`, `joint.fafmax.faf95_max_gen_anc`, and `joint.fafmax.faf95_max` from the gnomAD parquet table.

**Command:**
```bash
src/scripts/run_annotate_gnomad.sh output_clinvar.tsv output_final.tsv \
    --execution-mode athena \
    --gnomad-version v4.1 \
    --athena-output-location s3://your-bucket/athena-query-results/ \
    --athena-region us-east-1
```

**Environment variables (can be set in `settings/.env` instead of flags):**

| Variable | Flag | Description |
|---|---|---|
| `GNOMAD_ATHENA_DATABASE` | `--athena-database` | Athena database name (default: `gnomad`) |
| `GNOMAD_ATHENA_TABLE` | `--athena-table` | Table name (default: derived from `--gnomad-version`, e.g. `gnomad_joint_v4_1`) |
| `GNOMAD_ATHENA_OUTPUT_LOCATION` | `--athena-output-location` | S3 URI for Athena query result storage (**required**) |
| `GNOMAD_ATHENA_WORKGROUP` | `--athena-workgroup` | Optional Athena workgroup |
| `GNOMAD_ATHENA_REGION` or `AWS_REGION` | `--athena-region` | AWS region for the Athena client |
| `GNOMAD_ATHENA_ROW_BATCH_SIZE` | `--athena-row-batch-size` | Input rows processed per Athena query batch (default: 1000) |

**Additional tuning flags:**
- `--athena-max-caids-per-query N` — maximum number of CAIDs per Athena `IN` clause (default: 16250, matching MaveDB batch size)
- `--athena-poll-seconds N` — polling interval while waiting for Athena query completion (default: 5)

**Notes (Athena mode):**
- Requires `boto3` (included in project dependencies via `pyproject.toml`); AWS credentials must be available in the standard boto3 credential chain (environment variables, instance profile, `~/.aws/credentials`, etc.)
- Input rows are processed in batches; output is written and flushed after each batch, preserving input row order
- CAID lookups are cached in memory across batches to avoid redundant Athena queries within a single run
- Results are also persisted to Redis between runs when Redis is available (`GNOMAD_CACHE_REDIS_ENABLED`), avoiding Athena round-trips for variants seen in previous runs
- `--download-only` and `--refresh-cache` flags are ignored in Athena mode (no local cache involved)
- All six output columns are pipe-delimited and candidate-aligned (same format as Hail mode)

### Step 7: Annotate with SpliceAI Scores (Optional)

**Purpose:** Add SpliceAI splice-impact scores per DNA HGVS candidate.

See [docs/annotate_spliceai.md](docs/annotate_spliceai.md) for full reference documentation including HGVS parsing rules, indel normalization, and parallel lookup options.

**Input columns:** `mapped_hgvs_g` (from step 1 or step 2)

**Output columns:** `spliceai.ds_ag`, `spliceai.ds_al`, `spliceai.ds_dg`, `spliceai.ds_dl`, `spliceai.dp_ag`, `spliceai.dp_al`, `spliceai.dp_dg`, `spliceai.dp_dl`, `spliceai.max_delta_score`

**Precomputed mode (recommended):**
```bash
src/scripts/run_annotate_spliceai.sh output_clinvar.tsv output_spliceai.tsv \
  --mode precomputed \
  --precomputed-snv-vcf spliceai_scores.masked.snv.hg38.vcf.gz \
  --precomputed-indel-vcf spliceai_scores.masked.indel.hg38.vcf.gz
```

**One-time cache/index preparation (optional):**
```bash
src/scripts/run_annotate_spliceai.sh output_clinvar.tsv output_spliceai.tsv \
  --mode precomputed \
  --precomputed-snv-vcf spliceai_scores.masked.snv.hg38.vcf.gz \
  --precomputed-indel-vcf spliceai_scores.masked.indel.hg38.vcf.gz \
  --prepare-cache-only
```

**Compute mode (slow, requires local SpliceAI/TensorFlow install):**
```bash
src/scripts/run_annotate_spliceai.sh output_clinvar.tsv output_spliceai.tsv \
  --mode compute \
  --genome /usr/src/app/path/to/reference.fa.gz \
  --annotation grch38
```

**Notes:**
- In precomputed mode, source VCFs are copied into the SpliceAI cache volume and indexed (`.tbi`) if needed.
- For protein-derived rows with multiple DNA candidates, all SpliceAI output columns are pipe-delimited and position-aligned.
- `spliceai.max_delta_score` is `max(DS_AG, DS_AL, DS_DG, DS_DL)` for each candidate.
- Precomputed files can have coverage limitations (for example, missing classes of indels).

### Step 8: Annotate with ClinGen Expert Panel Classifications (Optional)

**Purpose:** Join each variant candidate against ClinGen Evidence Repository (erepo) expert-panel classifications using HGVS, ClinVar Variation ID, and/or CAID.

See [docs/annotate_erepo.md](docs/annotate_erepo.md) for full reference documentation including join-key details, cache management, and warning diagnostics.

**Input columns:** `mapped_hgvs_c` (from step 1 or step 2), optionally `dna_clingen_allele_id` and a ClinVar variation ID column

**Output columns:** `clingen_evidence_repository.Assertion`, `.Expert Panel`, `.Disease`, `.Mondo Id`, `.Mode of Inheritance`, `.Applied Evidence Codes (Met)`, `.Applied Evidence Codes (Not Met)`, `.Summary of interpretation`, `.ClinVar Variation Id`, `.Allele Registry Id`, `.PubMed Articles`, `.Guideline`, `.Approval Date`, `.Published Date`, `.Retracted`, `.Evidence Repo Link`, `.Uuid`, `.warnings`

**Command:**
```bash
src/scripts/run_annotate_erepo.sh input.tsv output.tsv
```

**Notes:**
- The full erepo classification TSV is downloaded once and cached locally; use `--refresh-cache` to re-download.
- Lookups use up to three join keys (HGVS, ClinVar Variation ID, CAID) by default; results are deduplicated by erepo UUID.
- For rows with multiple DNA candidates, all output columns are pipe-delimited and candidate-aligned.
- Only variants with a formal ClinGen expert-panel classification will be annotated; most variants will have empty output columns.

### Step 9: Annotate with VEP Mutational Consequence (Optional)

**Purpose:** Add a mutational consequence term from Ensembl VEP for each DNA variant row.

See [docs/annotate_vep.md](docs/annotate_vep.md) for full reference documentation including transcript selection logic, API call strategy, and Redis cache configuration.

**Input columns:** `mapped_hgvs_c` (preferred), `mapped_hgvs_g`, `mapped_hgvs_p`

**Output columns:** `vep.mutational_consequences`, `vep.most_severe_mutational_consequence`, `vep.consequence_source`, `vep.access_date`, `vep.error`

**Command:**
```bash
src/scripts/run_annotate_vep.sh annotated.tsv annotated_vep.tsv
```

**Notes:**
- Uses Ensembl REST API endpoints `/vep/human/hgvs` and `/variant_recoder/human`.
- For multi-candidate rows, each candidate is resolved independently; all output columns are pipe-delimited per candidate position.
- `vep.mutational_consequences` contains `^`-delimited consequence terms for transcript HGVS candidates, or the single most-severe term for genomic/protein inputs.
- `vep.most_severe_mutational_consequence` always contains a single term per candidate.
- For transcript `c.delins` candidates where transcript reference sequence equals ALT (`ref == alt`), output is normalized to `no_change` (including `vep.consequence_source=no_change`) and `vep.error` is cleared.
- `vep.error` is pipe-delimited per candidate and contains the API error message (e.g. `api_error:VEP HTTP 503: ...`) when a candidate's request failed; empty string otherwise.
- Output rows are streamed in input order and support `--skip` / `--limit`.

**Transcript selection:**

Ensembl VEP always resolves an input HGVS to a genomic coordinate and then annotates *every* transcript overlapping that position. The top-level `most_severe_consequence` field therefore reflects the worst outcome across all transcripts at the locus, not necessarily the consequence on the transcript referenced in the input HGVS string.

This script selects the consequence for the specific input transcript when possible:

| Input HGVS type | API flag | How the consequence is selected | `vep.consequence_source` |
|---|---|---|---|
| RefSeq transcript (`NM_`, `NR_`) | `refseq=1` | Matched by `transcript_id` in `transcript_consequences` | `transcript` |
| Ensembl transcript (`ENST`) | _(default)_ | Matched by `transcript_id` in `transcript_consequences` | `transcript` |
| Genomic (`NC_`, `chr:g.`) | _(default)_ | `most_severe_consequence` used — no specific transcript | `most_severe` |
| Protein (`NP_`) / recoder fallback | _(default)_ | `most_severe_consequence` used | `most_severe` |

When no matching transcript entry is found (e.g. the transcript version is not in Ensembl's RefSeq set), the script falls back to `most_severe_consequence` and sets `vep.consequence_source` to `most_severe`.

**Column priority matters:** Because transcript HGVS strings (`mapped_hgvs_c`) yield transcript-specific consequences, `mapped_hgvs_c` should appear before `mapped_hgvs_g` in `--hgvs-cols`. The default ordering `mapped_hgvs_c,mapped_hgvs_g,mapped_hgvs_p` reflects this.

**Redis caching:**

VEP API responses are cached in Redis as `(most_severe, all_consequences, source)` triples per HGVS string. On subsequent runs, cached results are used directly without re-querying Ensembl. Caching is enabled automatically when the Redis service is reachable and fails gracefully when it is not.

| Variable | Default | Description |
|---|---|---|
| `VEP_CACHE_ENABLED` | `true` | Set to `0` / `false` to disable |
| `VEP_CACHE_REDIS_URL` | `redis://redis:6379/0` | Redis connection URL (also falls back to `REDIS_URL`) |
| `VEP_CACHE_PREFIX` | `vep:v1` | Key namespace; bump the version suffix to invalidate the cache after a significant Ensembl release |
| `VEP_CACHE_TTL_SECONDS` | `86400` (1 day) | TTL for hits |
| `VEP_CACHE_MISS_TTL_SECONDS` | `86400` (1 day) | TTL for misses (no consequence returned) |

### Step 10: Annotate with MaveDB Functional Classifications (Optional)

**Purpose:** Classify each variant against the primary and investigator-provided calibrations for its MaveDB score set.

See [docs/annotate_mavedb.md](docs/annotate_mavedb.md) for full reference documentation including calibration selection logic, range-based vs. class-based classification, and caching behaviour.

**Input columns:** `variant_urn` (MaveDB variant URN), `score` (numeric variant score)

**Output columns:** `mavedb.primary_calibration.urn`, `.name`, `.url`, `.functional_class_label`, `.functional_classification`, `mavedb.investigator_provided_calibration.urn`, `.name`, `.url`, `.functional_class_label`, `.functional_classification`

**Command:**
```bash
src/scripts/run_annotate_mavedb.sh input.tsv output.tsv
```

**Notes:**
- The score set URN is extracted from the variant URN (the portion before `#`).
- Calibration lists are fetched once per unique score set and cached in memory for the run.
- When the primary and investigator-provided calibrations are the same object, both column groups hold identical values.
- For range-based calibrations, classification uses the numeric `score` column directly. For class-based calibrations, the variant is looked up by URN via the MaveDB API.
- Rows with no variant URN or no applicable calibration receive empty annotation columns.

### Step 11: Annotate with In-Silico Predictor Scores (Optional)

**Purpose:** Annotate variants with pre-computed missense pathogenicity scores from REVEL, AlphaMissense, and/or MutPred2. All tools operate on GRCh38 genomic positions and score missense SNVs only.

See [docs/annotate_predictors.md](docs/annotate_predictors.md) for full reference documentation including data file preparation, pipe-delimited column behaviour, and troubleshooting.

**Input columns:** `mapped_hgvs_g` (pipe-delimited genomic HGVS)

**Output columns:** `revel.score`, `alphamissense.pathogenicity`, `alphamissense.class`, `mutpred2.score` (only columns for supplied tools are written)

**Command:**
```bash
src/scripts/run_annotate_predictors.sh input.tsv output.tsv \
  --alphamissense-file AlphaMissense_hg38.tsv.gz
```

**Notes:**
- At least one predictor source must be provided: `--revel-file` or `--revel-cache-file` for REVEL, `--alphamissense-file` or `--alphamissense-cache-file` for AlphaMissense, or `--dbnsfp-file` for MutPred2.
- All currently supported predictors score missense SNVs only; other variant types receive empty annotation columns.
- REVEL and AlphaMissense scores are pipe-aligned to the input candidates. MutPred2 emits a single maximum score (protein-level model).
- Requires `tabix` (htslib) on `$PATH` only when a tabix-indexed file (`--revel-file`, `--alphamissense-file`, `--dbnsfp-file`) is configured.

### Step 12: Flatten DNA Variants (Optional)

**Purpose:** For annotated variant files with multi-candidate DNA variants (from reverse translation), produce a flattened output where each DNA candidate has its own row.

See [docs/flatten_dna_variants.md](docs/flatten_dna_variants.md) for full reference documentation including column auto-detection logic and troubleshooting. This is useful when you want one row per DNA variant instead of pipe-delimited lists.

**Input columns:** `mapped_hgvs_g`, `mapped_hgvs_c`, `dna_clingen_allele_id`, and any annotation columns (`spliceai.*`, `clinvar.*`, `gnomad.*`, `vep.*`)

**Output columns:** All input columns, with pipe-delimited columns expanded to one row per candidate

**Command:**
```bash
src/scripts/run_flatten_dna_variants.sh annotated.tsv dna_variants.tsv
```

**Explicit column specification (if auto-detection doesn't work):**
```bash
src/scripts/run_flatten_dna_variants.sh annotated.tsv dna_variants.tsv \
  --dna-variant-columns mapped_hgvs_g,mapped_hgvs_c,dna_clingen_allele_id
```

**Notes:**
- Drops all protein-only rows without DNA reverse translations (rows where all DNA variant columns are empty)
- Expands pipe-delimited columns: `mapped_hgvs_g`, `mapped_hgvs_c`, `dna_clingen_allele_id`, and any annotation columns
- Non-list columns (e.g., gene_symbol, raw_hgvs_pro) are repeated across expanded rows
- Useful for downstream analysis that requires one genomic variant per row
- If all rows are protein-only without DNA variants, the script will error and produce no output

**Example input (with multi-candidate row):**
```
raw_hgvs_nt    raw_hgvs_pro    mapped_hgvs_g              mapped_hgvs_c                  dna_clingen_allele_id    spliceai.max_delta_score
              p.Arg175His     g.1A>T|g.2A>T|g.3A>T      c.1A>T|c.2A>T|c.3A>T          CA1|CA2|CA3             0.5|0.7|0.3
```

**Example output (3 rows):**
```
raw_hgvs_nt    raw_hgvs_pro    mapped_hgvs_g    mapped_hgvs_c    dna_clingen_allele_id    spliceai.max_delta_score
              p.Arg175His     g.1A>T           c.1A>T           CA1                      0.5
              p.Arg175His     g.2A>T           c.2A>T           CA2                      0.7
              p.Arg175His     g.3A>T           c.3A>T           CA3                      0.3
```


### Complete Pipeline Example

Processing a file end-to-end:

```bash
# Start with raw variants
src/scripts/run_map_variants.sh variants_raw.tsv variants_mapped.tsv

# Reverse-translate protein variants into DNA candidates
src/scripts/run_reverse_translate_protein_variants.sh variants_mapped.tsv variants_rt.tsv

# Add DNA-level ClinGen allele IDs
src/scripts/run_add_dna_clingen_allele_ids.sh variants_rt.tsv variants_clingen.tsv

# Add parsed position/allele columns (optional)
src/scripts/run_add_vcf_identifiers.sh variants_clingen.tsv variants_parsed.tsv

# Annotate with ClinVar (optional)
src/scripts/run_annotate_clinvar.sh variants_parsed.tsv variants_clinvar.tsv

# Annotate with gnomAD (optional, first time includes download)
src/scripts/run_annotate_gnomad.sh variants_clinvar.tsv variants_final.tsv

# Annotate with SpliceAI (optional)
src/scripts/run_annotate_spliceai.sh variants_final.tsv variants_spliceai.tsv \
  --mode precomputed \
  --precomputed-snv-vcf spliceai_scores.masked.snv.hg38.vcf.gz \
  --precomputed-indel-vcf spliceai_scores.masked.indel.hg38.vcf.gz

# Annotate with VEP mutational consequence (optional)
src/scripts/run_annotate_vep.sh variants_spliceai.tsv variants_vep.tsv

# Flatten multi-candidate variants to one row per DNA variant (optional)
src/scripts/run_flatten_dna_variants.sh variants_vep.tsv variants_dna_only.tsv
```

### Data Volume Management

ClinVar, gnomAD, and SpliceAI annotation steps use Docker volumes to persist caches:

```yaml
# In compose.yaml

volumes:
```
  variant-annotation-clinvar-cache:
  variant-annotation-gnomad-cache:
  variant-annotation-spliceai-cache:

services:
  annotate-clinvar:
    volumes:
      - variant-annotation-clinvar-cache:/clinvar-cache

  annotate-gnomad:
    volumes:
      - variant-annotation-gnomad-cache:/gnomad-cache

  annotate-spliceai:
    volumes:
      - variant-annotation-spliceai-cache:/spliceai-cache

  annotate-vep:
    volumes:
      - ${VARIANT_DATA_DIR:-./data}:/work
```

This ensures caches survive container restarts and are shared across runs.

## Redis Caching

Several pipeline scripts use Redis to cache API responses between runs, avoiding redundant network requests when the same variant is processed multiple times (e.g. across incremental runs or overlapping datasets). All caching degrades gracefully: if Redis is unreachable or disabled, scripts continue without caching.

The Redis service is started automatically by Docker Compose and is available at `redis://redis:6379/0` inside the container network.

### ClinGen Allele Registry

**Used by:** `map_variants` (steps 1 & 3 combined), `add_dna_clingen_allele_ids`, `annotate_clinvar`
**Implementation:** `src/lib/clingen.py`

Two types of ClinGen requests are cached:

| Request type | Endpoint | Redis key | Value stored |
|---|---|---|---|
| HGVS lookup | `GET /allele?hgvs=<HGVS>` | `clingen:v1:hgvs:<HGVS>` | Allele ID string (e.g. `CA123456`) |
| Allele ID lookup | `GET /allele/<CA_ID>` | `clingen:v1:allele:<CA_ID>` | Full JSON response body |

A successful HGVS lookup writes **both** keys: the mapping key (`hgvs:`) pointing at the allele ID, and the allele key (`allele:`) holding the complete response. The full response is cached rather than just the allele ID because multiple callers extract different fields from the same record — `annotate_clinvar` needs `externalRecords.ClinVarVariations` and `externalRecords.ClinVarAlleles`, while coordinate-resolution helpers need `genomicAlleles`. Caching the full body under the allele ID means any caller, whether it arrived via HGVS or directly by allele ID, can be served from a single Redis entry.

The two-key design also handles partial eviction gracefully: if the `hgvs:` key is still present but the `allele:` key was evicted, the code detects this and re-fetches by allele ID without re-querying by HGVS.

`map_variants` additionally caches `dcd_mapping` HGVS normalization results:

| Redis key | Value stored |
|---|---|
| `clingen:v1:genomic_hgvs:<HGVS>` | Normalized genomic HGVS string |

This key is populated by a patch applied to `dcd_mapping`'s internal `fetch_clingen_genomic_hgvs()` call, which itself queries ClinGen during variant normalization for accession-referenced inputs.

Both hits and misses (404s, placeholder IDs) are cached. Misses are stored as the sentinel `__MISS__`. Default TTL is 86,400 s (1 day) for both.

| Environment variable | Default | Description |
|---|---|---|
| `CLINGEN_CACHE_ENABLED` | `true` | Set to `0`/`false` to disable |
| `CLINGEN_CACHE_REDIS_URL` | `redis://redis:6379/0` | Falls back to `REDIS_URL` |
| `CLINGEN_CACHE_PREFIX` | `clingen:v1` | Bump version suffix to invalidate all entries |
| `CLINGEN_CACHE_TTL_SECONDS` | `86400` | TTL for hit entries |
| `CLINGEN_CACHE_MISS_TTL_SECONDS` | `86400` | TTL for miss sentinels |

Use `src/scripts/run_clear_clingen_cache.sh` to delete all ClinGen cache keys and force fresh API lookups.

### VEP (Ensembl REST API)

**Used by:** `annotate_vep`
**Implementation:** `src/annotate_vep.py`

VEP consequences are cached per HGVS string. Each entry stores the full consequence result as a small JSON object:

| Redis key | Value stored |
|---|---|
| `vep:v1:<HGVS>` | `{"c": "<most_severe>", "cs": ["<term>", ...], "s": "<source>"}` |

Where `c` is the most-severe consequence, `cs` is the full list of consequences from the matched transcript (or `null` for genomic/protein inputs), and `s` is the source tag (`"transcript"`, `"most_severe"`, etc.). Storing the full consequence list avoids a second API call if downstream code needs it.

Results that failed with a transient API error (`source == "api_error"`) are **not** cached, so they are retried on the next run. All other results — including definitive misses — are cached with the sentinel `__MISS__`.

Reads and writes use Redis pipeline operations (`MGET`/`SET`) to minimise round-trips when processing large batches.

| Environment variable | Default | Description |
|---|---|---|
| `VEP_CACHE_ENABLED` | `true` | Set to `0`/`false` to disable |
| `VEP_CACHE_REDIS_URL` | `redis://redis:6379/0` | Falls back to `REDIS_URL` |
| `VEP_CACHE_PREFIX` | `vep:v1` | Bump version suffix after a major Ensembl release |
| `VEP_CACHE_TTL_SECONDS` | `86400` | TTL for hit entries |
| `VEP_CACHE_MISS_TTL_SECONDS` | `86400` | TTL for miss sentinels |

### gnomAD

**Used by:** `annotate_gnomad` (both Hail and Athena modes)
**Implementation:** `src/annotate_gnomad.py`

Full gnomAD records are cached per ClinGen allele ID (CAID), or per genomic coordinate string in coordinate-lookup mode:

| Redis key | Value stored |
|---|---|
| `gnomad:v1:<CAID>` | JSON-serialized `GnomadRecord` |

Each cached record includes allele counts and frequencies, filtering allele frequencies, per-callset QC filter flags, VEP gene symbols, and any histogram data (age/allele-balance distributions) that was present when the record was fetched. Storing the full record avoids re-querying the Hail table or Athena for any field combination, regardless of which `--age-histograms` / `--allele-balance-histograms` flags the current run uses.

Only successful lookups are cached; variants not found in gnomAD are not stored (no miss sentinel). The default TTL is 604,800 s (7 days), reflecting that gnomAD releases are infrequent.

Reads and writes use Redis pipeline operations for batch efficiency.

| Environment variable | Default | Description |
|---|---|---|
| `GNOMAD_CACHE_REDIS_ENABLED` | `true` | Set to `0`/`false` to disable |
| `GNOMAD_CACHE_REDIS_URL` | `redis://redis:6379/0` | Falls back to `REDIS_URL` |
| `GNOMAD_CACHE_REDIS_PREFIX` | `gnomad:v1` | Bump version suffix after a gnomAD release |
| `GNOMAD_CACHE_REDIS_TTL_SECONDS` | `604800` | TTL for all cached entries (7 days) |

### Summary

| Script | Key prefix | Cached value | TTL (hit/miss) |
|---|---|---|---|
| ClinGen (HGVS lookup) | `clingen:v1:hgvs:` | Allele ID string | 1 day / 1 day |
| ClinGen (allele record) | `clingen:v1:allele:` | Full JSON response | 1 day / 1 day |
| ClinGen (genomic HGVS norm.) | `clingen:v1:genomic_hgvs:` | Normalized HGVS string | 1 day / 1 day |
| VEP | `vep:v1:` | Consequence + source JSON | 1 day / 1 day |
| gnomAD | `gnomad:v1:` | Full GnomadRecord JSON | 7 days / not cached |

## Bundled Utilities

The following tools are included alongside the main pipeline for data inspection, cross-validation, and TSV manipulation. They are not annotation steps.

| Script | Description |
|---|---|
| `run_backfill_clingen_allele_id.sh` | Backfill missing `clingen_allele_id` values (Step 1 output only) where the initial `map_variants` run left cells blank. See [docs/backfill_clingen_allele_ids.md](docs/backfill_clingen_allele_ids.md) |
| `run_clear_clingen_cache.sh` | Delete ClinGen Allele Registry Redis cache keys to force fresh API lookups. See [docs/clear_clingen_cache.md](docs/clear_clingen_cache.md) |
| `run_compare_cvfg_datasets.sh` | Cross-validate a CVFG pipeline flat file against the integrated variant effect dataset via a positional-genome join; produces matched, unmatched, and error output files. See [docs/compare_cvfg_datasets.md](docs/compare_cvfg_datasets.md) |
| `run_diff_hgvs.sh` | Compare `mapped_hgvs_g/c/p` columns between two TSV versions, reporting each changed pipe-component. See [docs/diff_hgvs.md](docs/diff_hgvs.md) |
| `run_diff_vcf_coords.sh` | Compare VCF-style coordinate columns (genomic, transcript, protein) between two TSV versions. See [docs/diff_vcf_coords.md](docs/diff_vcf_coords.md) |
| `run_diff_tsv.sh` | Row-level diff between two TSV versions: changed rows and HGVS component diffs. See [docs/diff_tsv.md](docs/diff_tsv.md) |
| `run_haplotypes_to_delins.sh` | Rewrite intra-codon c.-haplotype HGVS expressions to delins notation. See [docs/haplotypes_to_delins.md](docs/haplotypes_to_delins.md) |
| `run_gnomad_hail_shell.sh` | Launch an interactive bash shell in the Hail-enabled gnomAD container for ad-hoc table inspection. See [docs/gnomad_hail_shell.md](docs/gnomad_hail_shell.md) |
| `run_filter_dna_variants.sh` | Remove non-SNV DNA candidates while preserving row count and aligned candidate columns. See [docs/filter_dna_variants.md](docs/filter_dna_variants.md) |
| `run_utilities.sh` | General-purpose TSV manipulation: `filter-rows`, `filter-columns`, `merge-columns`, `rename-columns`, `reorder-columns`, `replace-rows`, `compare-columns`. See [docs/utilities.md](docs/utilities.md) |

## Pipeline Data Flow Diagram

```
                        INPUT (raw variants)
                               ↓
                    ┌─────────────────────┐
                    │   map_variants      │ ✓ Required
                    │                     │
                    │ Output:             │
                    │ mapped_hgvs_g/c/p   │
                    │ clingen_allele_id   │
                    └─────────────────────┘
                               ↓
              ┌───────────────────────────────┐
              │ Is variant protein-only?      │
              │ (no c./g., yes p.)            │
              └────┬─────────────────┬────────┘
                   │ YES             │ NO
                   ↓                 │
        ┌──────────────────────┐     │
        │ reverse_translate    │     │
        │ _protein_variants    │     │ (DNA variants pass through)
        │                      │     │
        │ Output: candidates   │     │
        │ in c./g. format      │     │
        └──────────────────────┘     │
                   │                 │
                   └─────────────┬───┘
                                 ↓
                   ┌─────────────────────────────┐
                   │ add_dna_clingen_allele_ids  │ ✓ Required for
                   │                             │   annotation
                   │ Output:                     │
                   │ dna_clingen_allele_id       │
                   │ (pipe-delimited candidates) │
                   └─────────────────────────────┘
                                 ↓
              ┌──────────────────────────────────┐
              │ Want parsed positions/alleles?   │
              └────┬─────────────────┬───────────┘
                   │ YES             │ NO
                   ↓                 │
        ┌────────────────────────┐   │
        │ add_vcf_identifiers    │   │
        └────────────────────────┘   │
                   │                 │
                   └─────────────┬───┘
                                 ↓
              ┌──────────────────────────────────┐
              │ Want ClinVar annotations?        │
              └────┬─────────────────┬───────────┘
                   │ YES             │ NO
                   ↓                 │
        ┌─────────────────────────┐  │
        │ annotate_clinvar        │  │
        │                         │  │
        │ Caches: NCBI ClinVar    │  │
        │ TSV (monthly)           │  │
        └─────────────────────────┘  │
                   │                 │
                   └─────────────┬───┘
                                 ↓
              ┌──────────────────────────────────┐
              │ Want gnomAD allele frequencies?  │
              └────┬─────────────────┬───────────┘
                   │ YES             │ NO
                   ↓                 │
        ┌─────────────────────────┐  │
        │ annotate_gnomad         │  │
        │                         │  │
        │ Caches: Hail table      │  │
        │ (indexed, keyed by CAID)│  │
        └─────────────────────────┘  │
                   │                 │
                   └─────────────┬───┘
                                 ↓
              ┌──────────────────────────────────┐
              │ Want SpliceAI splice scores?     │
              └────┬─────────────────┬───────────┘
                   │ YES             │ NO
                   ↓                 │
        ┌─────────────────────────┐  │
        │ annotate_spliceai       │  │
        │                         │  │
        │ Caches: precomputed     │  │
        │ VCFs or uses compute    │  │
        │ mode with TensorFlow    │  │
        └─────────────────────────┘  │
                   │                 │
                   └─────────────┬───┘
                                 ↓
              ┌──────────────────────────────────┐
              │ Want one row per DNA variant?    │
              │ (flatten multi-candidate rows)   │
              └────┬─────────────────┬───────────┘
                   │ YES             │ NO
                   ↓                 │
        ┌─────────────────────────┐  │
        │ flatten_dna_variants    │  │
        │                         │  │
        │ Drops protein-only rows │  │
        │ (no DNA variants)       │  │
        └─────────────────────────┘  │
                   │                 │
                   └─────────────┬───┘
                                 ↓
                       OUTPUT (fully annotated variants)
```

## Troubleshooting Pipeline Issues

### Step 1: map_variants Issues

**Problem: Variants not mapping to genomic coordinates**
- Check that `raw_hgvs_nt` or `raw_hgvs_pro` columns are present and non-empty
- For sequence-based mapping, ensure `target_sequence` is present
- For transcript-based mapping (e.g., `NM_`:), verify the transcript exists in the reference genome
- Check the `mapped_errors` column in output for specific error details

**Problem: BLAT fails with error 137 (Out of Memory)**

The BLAT subprocess was killed by the operating system due to out-of-memory (OOM) conditions. This typically occurs when processing very large groups of variants in a single BLAT run, as each variant pair requires significant memory for alignment.

The `map-variants` command automatically enables error-driven retry with progressive chunking by default. If error 137 still persists:

1. **Reduce initial chunk size** (default: 500):
   ```bash
   src/scripts/run_map_variants.sh input.tsv output.tsv --dcd-chunk-size-on-137 250
   ```

2. **Increase maximum retry attempts** (default: 3):
   ```bash
   src/scripts/run_map_variants.sh input.tsv output.tsv --dcd-max-retry-attempts 5
   ```

3. **Disable retry and use manual chunking** (only if needed):
   ```bash
   src/scripts/run_map_variants.sh input.tsv output.tsv --no-dcd-chunk-on-137
   ```
   Then pre-split your input file by gene/target sequence group before rerunning.

4. **Increase Docker memory limit**:
  Edit `compose.yaml` and add memory limits under the `map-variants` service:
   ```yaml
   services:
     map-variants:
       mem_limit: 8g  # or higher as needed
       memswap_limit: 8g
   ```

5. **Create a retry input file** with only failed rows from a prior run:
   - The output file will contain error rows with the mapping error column populated
   - Filter to only error 137 rows: `error_code_137.tsv`
   - Rerun on this smaller file with aggressive chunking:
   ```bash
   src/scripts/run_map_variants.sh error_code_137.tsv output_retry.tsv \
       --dcd-chunk-size-on-137 100 \
       --dcd-max-retry-attempts 5
   ```

For more details on chunking options, see the [BLAT Error 137 Retry Strategy](#blat-error-137-retry-strategy) section earlier in this document.

**Problem: Slow mapping performance**
- Use `--merge-existing` to reuse prior results (avoid remapping identical variants)
- Filter input to only unmapped variants before running
- Consider splitting large files by gene/target sequence with `--group-by`

### Step 2: reverse_translate_protein_variants Issues

**Problem: Reverse translation produces no candidates**
- Ensure `mapped_hgvs_p` (protein HGVS) is present and well-formed (e.g., `p.Ala406Thr`)
- Check that `mapped_hgvs_c` and `mapped_hgvs_g` are empty (DNA variants won't be re-translated)
- Verify the UTA database is running: `docker compose logs uta`
- Check for codon degeneracy; some protein changes have very few DNA backtranslations

**Problem: Reverse translation is slow**
- This is expected; UTA lookups can be time-consuming on large sets
- Consider filtering to only protein-only rows before running

### Step 3: add_dna_clingen_allele_ids Issues

**Problem: ClinGen lookup fails (many empty results)**
- ClinGen API may be rate-limited; script includes automatic retry with exponential backoff
- Check internet connectivity to `reg.genome.network`
- Verify HGVS strings are well-formed (check prior steps' output)
- Some valid variants may not exist in ClinGen registry

**Problem: Inconsistent pipe-delimited output lengths**
- This is normal; empty slots (e.g., `CA1||CA3`) indicate candidates without ClinGen matches
- Check earlier reverse-translation step for the candidate count per row

### Step 4: add_vcf_identifiers Issues

**Problem: Parsing errors or missing ref alleles**
- UTA database may not be running; start it: `docker compose up uta`
- Some variants may have incomplete HGVS strings; safe to skip these rows
- Check logs for specific HGVS parsing errors

### Step 5: annotate_clinvar Issues

**Problem: ClinVar lookup returns no annotations (empty columns)**
- First-run downloads the NCBI ClinVar TSV; check cache directory is writable
- ClinGen → ClinVar ID lookup may fail for valid ClinGen alleles not in ClinVar
- Try using a different `--clinvar-version` (e.g., older monthly archive)

**Problem: ClinVar cache is stale or corrupt**
- Delete the cache directory: `rm -rf ./clinvar_cache`
- Re-run; a fresh cache will be downloaded
- Specify `--cache-dir` explicitly if default location doesn't work

**Problem: Docker permission denied on cache volume**
- Ensure `variant-annotation-clinvar-cache` volume exists: `docker volume ls | grep clinvar`
- If missing, start Docker Compose once: `docker compose up -d`

### Step 6: annotate_gnomad Issues

**Problem: "Hail is required for gnomAD annotation"**
- gnomAD extra not installed; reinstall: `pip install -e '.[gnomad]'`
- Ensure Java is available on the system (required by Hail)
- In Docker, restart and rebuild: `docker compose up --build`

**Problem: Hail table download stalls or fails**
- First-time download from GCP public data can take 10–30 minutes
- Check internet connectivity and GCP access
- Use `--download-only` flag to debug download separately: `src/scripts/run_annotate_gnomad.sh dummy.tsv dummy_out.tsv --download-only`
- If download fails, use `--refresh-cache` on next attempt

**Problem: gnomAD lookup returns no matches despite valid CAIDs**
- Not all ClinGen alleles exist in gnomAD; some are too rare
- Check that `dna_clingen_allele_id` is correctly populated from step 3
- Try a different gnomAD version: `--gnomad-version v4.0`

**Problem: Docker volume mount issues with gnomAD cache**
- Verify volume exists: `docker volume ls | grep gnomad`
- Ensure volume has write permissions
- Check Docker disk space; Hail table cache can be several GB

## Performance and Optimization Tips

### General Strategies

**Batch processing**
- Process variants in batches of 10,000–50,000 rows at a time
- Smaller batches allow easier error recovery if a step fails mid-run
- Larger batches are more efficient for network requests (ClinGen, NCBI, GCP)

**Reuse cached results**
- For `map_variants`: use `--merge-existing` to skip remapped variants
- For `annotate_clinvar`: cache persists in Docker volume across runs
- For `annotate_gnomad`: first run is slow; subsequent runs are ~10x faster (local cache hits)

**Monitor resource usage**
- Watch Docker memory: `docker stats`
- On Apple Silicon, `map-variants` uses emulation; slower than native runs
- Increase Docker CPU/memory allocation in settings if hitting limits

### Step-Specific Optimizations

**map_variants**
- Group by gene or target sequence to align BLAT calls
- Merge prior results to avoid reprocessing: `--merge-existing output_v1.tsv --merge-match-col gene_symbol`
- For very large files (>100k variants), pre-filter to only unmapped rows

**reverse_translate_protein_variants**
- Only runs on protein-only rows; DNA variants skip this step
- Consider filtering input to only protein-only rows to avoid unnecessary I/O
- UTA lookups are I/O-bound; parallelization isn't feasible within a single run

**add_dna_clingen_allele_ids**
- ClinGen API rate limits to ~1 req/sec; handle gracefully (script includes backoff)
- Caches HGVS lookups in-memory to avoid duplicate requests
- For repeated runs on similar data, consider externalizing this cache

**annotate_clinvar**
- NCBI ClinVar TSV is ~20–50 MB; download is fast (< 1 min)
- Lookup is O(1) hash after parsing; very fast for large result sets
- Docker volume ensures cache survives restarts

**annotate_gnomad**
- First-time Hail table download: 10–30 minutes (network-dependent)
- Second+ runs: <5 seconds (local cache indexed read)
- Hail indexing: ~5 minutes (one-time per table version)
- Cache can be several GB; ensure sufficient disk space

## Example Input and Output Data

### Example 1: DNA-only Variant (Simple Case)

**Input to Step 1 (map_variants):**
```
gene_symbol	raw_hgvs_nt	raw_hgvs_pro	target_sequence
BRCA1	NM_007294.3:c.68_69delAG	p.Glu23fs	
```

**Output of Step 1:**
```
gene_symbol	raw_hgvs_nt	raw_hgvs_pro	target_sequence	mapped_hgvs_g	mapped_hgvs_c	mapped_hgvs_p	clingen_allele_id	mapped_errors
BRCA1	NM_007294.3:c.68_69delAG	p.Glu23fs		g.43044394_43044395delAG	NM_007294.3:c.68_69delAG	p.Glu23fs	CA324372	
```

**Output of Step 3 (add_dna_clingen_allele_ids):**
```
gene_symbol	raw_hgvs_nt	raw_hgvs_pro	mapped_hgvs_g	mapped_hgvs_c	mapped_hgvs_p	clingen_allele_id	dna_clingen_allele_id
BRCA1	NM_007294.3:c.68_69delAG	p.Glu23fs	g.43044394_43044395delAG	NM_007294.3:c.68_69delAG	p.Glu23fs	CA324372	CA324372
```

### Example 2: Protein-Only Variant with Multiple Candidates (Complex Case)

**Input to Step 1:**
```
gene_symbol	raw_hgvs_nt	raw_hgvs_pro	target_sequence
TP53		p.Arg175His	MDDLMLSPDD...
```

**Output of Step 1:**
```
gene_symbol	raw_hgvs_nt	raw_hgvs_pro	target_sequence	mapped_hgvs_g	mapped_hgvs_c	mapped_hgvs_p	clingen_allele_id	mapped_errors
TP53		p.Arg175His	MDDLMLSPDD...		p.Arg175His			
```

**After Step 2 (reverse_translate_protein_variants):**
```
gene_symbol	raw_hgvs_nt	raw_hgvs_pro	mapped_hgvs_g	mapped_hgvs_c	mapped_hgvs_p	mapped_errors
TP53		p.Arg175His	g.7571720C>G|g.7571720C>A|g.7571721G>A	c.524C>G|c.524C>A|c.525G>A	p.Arg175His	
```
(Protein Arg175His maps to 3 different DNA changes due to codon degeneracy)

**After Step 3 (add_dna_clingen_allele_ids):**
```
gene_symbol	mapped_hgvs_g	mapped_hgvs_c	mapped_hgvs_p	dna_clingen_allele_id
TP53	g.7571720C>G|g.7571720C>A|g.7571721G>A	c.524C>G|c.524C>A|c.525G>A	p.Arg175His	CA123456|CA123457||
```
(Note: 3rd candidate has no ClinGen match, shown as empty slot)

**After Step 5 (annotate_clinvar):**
```
gene_symbol	dna_clingen_allele_id	clinvar.202601.clinical_significance	clinvar.202601.review_status	clinvar.202601.stars	clinvar.202601.last_review_date
TP53	CA123456|CA123457||	Pathogenic	reviewed by expert panel	3	2023-06-15
```
(Uses first successful candidate: CA123456 found in ClinVar)

**After Step 6 (annotate_gnomad):**
```
gene_symbol	dna_clingen_allele_id	clinvar.202601.clinical_significance	gnomad.v4.1.allele_frequency	gnomad.v4.1.minor_allele_frequency	gnomad.v4.1.allele_count
TP53	CA123456|CA123457||	Pathogenic	0.00234	0.00234	1547
```

## Complete Column Reference Guide

### Input Columns (User-Provided)

| Column | Type | Required | Description |
|--------|------|----------|-------------|
| `raw_hgvs_nt` | string | Conditional | Nucleotide HGVS variant (c. or n. format, optionally with transcript prefix like `NM_:`) |
| `raw_hgvs_pro` | string | Conditional | Protein HGVS variant (p. format, e.g., `p.Ala406Thr`) |
| `target_sequence` | string | Conditional | Nucleotide or amino acid sequence; required for sequence-based mapping |
| `<group-by-col>` | string | Optional | Column to group by during mapping (e.g., `gene_symbol`, `target_sequence` hash) |

**Conditional Requirements:**
- At least one of `raw_hgvs_nt` or `raw_hgvs_pro` must be present per row
- For sequence-based mapping (no transcript prefix), `target_sequence` is required

---

### Output Columns (Generated by Pipeline)

#### Step 1: map_variants

| Column | Type | Description |
|--------|------|-------------|
| `mapped_hgvs_g` | string | Genomic HGVS on GRCh38 (g. format) |
| `mapped_hgvs_c` | string | Transcript HGVS (c. format, e.g., `NM_007294.3:c.68_69delAG`) |
| `mapped_hgvs_p` | string | Protein HGVS (p. format) |
| `clingen_allele_id` | string | ClinGen Allele ID from Registry (e.g., `CA123456`) |
| `mapped_errors` | string | Error message if mapping failed (e.g., "Transcript not found") |

#### Step 2: reverse_translate_protein_variants

| Column | Type | Description |
|--------|------|-------------|
| `mapped_hgvs_c` | string | Updated with pipe-delimited c. candidates (e.g., `c.1218G>A\|c.1221G>A`) or unchanged if DNA |
| `mapped_hgvs_g` | string | Updated with pipe-delimited g. candidates or unchanged if DNA |

#### Step 3: add_dna_clingen_allele_ids

| Column | Type | Description |
|--------|------|-------------|
| `dna_clingen_allele_id` | string | Pipe-delimited ClinGen Allele IDs (e.g., `CA123456\|\|CA123457`); empty slots for unresolved candidates |

#### Step 4: add_vcf_identifiers

| Column | Type | Description |
|--------|------|-------------|
| `mapped_hgvs_g_start` | int | Start position of genomic variant |
| `mapped_hgvs_g_stop` | int | End position of genomic variant |
| `mapped_hgvs_g_ref` | string | Reference allele (genomic) |
| `mapped_hgvs_g_alt` | string | Alternate allele (genomic) |
| `mapped_hgvs_c_start` | int | Start position of transcript variant |
| `mapped_hgvs_c_stop` | int | End position of transcript variant |
| `mapped_hgvs_c_ref` | string | Reference allele (transcript) |
| `mapped_hgvs_c_alt` | string | Alternate allele (transcript) |
| `mapped_hgvs_p_start` | int | Start position of protein variant |
| `mapped_hgvs_p_stop` | int | End position of protein variant |
| `mapped_hgvs_p_ref` | string | Reference amino acid |
| `mapped_hgvs_p_alt` | string | Alternate amino acid |
| `touches_intronic_region` | string | Pipe-delimited per-candidate flags (`"true"`/`"false"`), aligned to `mapped_hgvs_c` |
| `spans_intron` | string | Pipe-delimited per-candidate flags (`"true"`/`"false"`), aligned to `mapped_hgvs_c` |

#### Step 5: annotate_clinvar

| Column | Type | Description |
|--------|------|-------------|
| `clinvar.<YYYYMM>.clinical_significance` | string | ClinVar significance (e.g., `Pathogenic`, `Benign`, `Uncertain significance`) |
| `clinvar.<YYYYMM>.review_status` | string | Review status (e.g., `reviewed by expert panel`, `criteria provided`) |
| `clinvar.<YYYYMM>.stars` | int | Star rating: 0–4 (higher = more confident) |
| `clinvar.<YYYYMM>.last_review_date` | string | ISO 8601 date (e.g., `2023-06-15`) |

(Default namespace is `clinvar`; customize with `--namespace` flag)

#### Step 6: annotate_gnomad

Standard columns (always emitted):

| Column | Type | Description |
|--------|------|-------------|
| `gnomad.<VERSION>.minor_allele_frequency` | float | MAF = min(AF, 1-AF) |
| `gnomad.<VERSION>.allele_frequency` | float | Allele frequency in gnomAD cohort |
| `gnomad.<VERSION>.allele_count` | int | Count of alternate alleles in cohort |
| `gnomad.<VERSION>.allele_number` | int | Total alleles in cohort (2× sample count) |
| `gnomad.<VERSION>.faf95_max` | float | Filtering allele frequency (95% CI max) |
| `gnomad.<VERSION>.faf95_max_ancestry` | string | Ancestry group with highest FAF95 |
| `gnomad.<VERSION>.filters` | string | Pipe-delimited combined QC filter flags (empty = PASS) |
| `gnomad.<VERSION>.exome_filters` | string | Pipe-delimited exome-callset QC filter flags (empty = PASS or not in exome callset) |
| `gnomad.<VERSION>.genome_filters` | string | Pipe-delimited genome-callset QC filter flags (empty = PASS or not in genome callset) |
| `gnomad.<VERSION>.gene_symbols` | string | Pipe-delimited VEP gene symbols overlapping the variant |

Optional histogram columns (Hail mode only, enabled by `--age-histograms` / `--allele-balance-histograms`):

| Column pattern | Type | Description |
|--------|------|-------------|
| `gnomad.<V>.<hist_name>.n_smaller` | int | Count of carriers below the lowest bin edge |
| `gnomad.<V>.<hist_name>.bin_N_<lo>_<hi>` | int | Carrier count in bin N with edges \[lo, hi) |
| `gnomad.<V>.<hist_name>.n_larger` | int | Count of carriers above the highest bin edge |

Where `<hist_name>` is one of: `age_hist_exome_het`, `age_hist_exome_hom`, `age_hist_genome_het`, `age_hist_genome_hom`, `age_hist_joint_het`, `age_hist_joint_hom` (age histograms) or `ab_hist_exome_adj`, `ab_hist_exome_raw`, `ab_hist_genome_adj`, `ab_hist_genome_raw`, `ab_hist_joint_adj`, `ab_hist_joint_raw` (allele balance histograms). See [docs/annotate_gnomad.md](docs/annotate_gnomad.md) for full details.

(Default version is `v4.1`; customize with `--gnomad-version` flag)

#### Step 7: annotate_spliceai

| Column | Type | Description |
|--------|------|-------------|
| `spliceai.ds_ag` | float | Acceptor gain delta score |
| `spliceai.ds_al` | float | Acceptor loss delta score |
| `spliceai.ds_dg` | float | Donor gain delta score |
| `spliceai.ds_dl` | float | Donor loss delta score |
| `spliceai.dp_ag` | float | Acceptor gain delta position |
| `spliceai.dp_al` | float | Acceptor loss delta position |
| `spliceai.dp_dg` | float | Donor gain delta position |
| `spliceai.dp_dl` | float | Donor loss delta position |
| `spliceai.max_delta_score` | float | Max of DS_AG, DS_AL, DS_DG, DS_DL |

#### Step 11: flatten_dna_variants

**No new columns are added in Step 11.** Instead, pipe-delimited columns from previous steps are expanded so that each DNA candidate gets its own row. The output file contains all columns from the input file, but with:
- One row per DNA candidate (instead of one row per protein with pipe-delimited candidates)
- Non-list columns (e.g., `raw_hgvs_pro`, `gene_symbol`) repeated across the expanded rows
- Protein-only rows (without DNA variants) dropped entirely

This step is a format transformation, not an annotation enrichment.

---

### Key Properties of Pipe-Delimited Columns

When a row has multiple DNA candidates (from reverse translation), the following columns are pipe-delimited:

**Core DNA variant columns:**
- `mapped_hgvs_c`
- `mapped_hgvs_g`
- `dna_clingen_allele_id`

**Parsed position/allele columns (from Step 4):**
- `mapped_hgvs_g_start`, `mapped_hgvs_g_stop`, `mapped_hgvs_g_ref`, `mapped_hgvs_g_alt`
- `mapped_hgvs_c_start`, `mapped_hgvs_c_stop`, `mapped_hgvs_c_ref`, `mapped_hgvs_c_alt`
- `touches_intronic_region`, `spans_intron`

**Annotation columns:**
- `spliceai.*` output columns (after step 7)
- `clinvar.*` output columns (after step 5)
- `gnomad.*` output columns (after step 6)

**Important:** Empty slots are preserved (e.g., `CA1||CA3` has 3 candidates, 2nd without a match).
This ensures positional alignment across all downstream annotations.

When annotating, the pipeline tries candidates in order and returns the first successful match.

**Note on Step 11 (flatten_dna_variants):** This optional step converts pipe-delimited columns into separate rows, with one row per DNA candidate. After flattening, there are no more pipe-delimited values in these columns—each candidate has its own row.

## Local installation (without Docker)

Install the package in editable mode with all extras:

```bash
pip install -e '.[dev,tests]'
```

Install dcd_mapping from GitHub without dependencies (to avoid a transient ga4gh.vrs pin conflict):

```bash
pip install --no-deps 'dcd-mapping @ git+https://github.com/VariantEffect/dcd_mapping2.git@main'
```

Install only the runtime dependencies:

```bash
pip install -e .
```


