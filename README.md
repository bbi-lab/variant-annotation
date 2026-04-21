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

This step is required for transcript lookups in `reverse_translate_protein_variants` and `add_variant_position_alleles`. It downloads the dump to a volume-mounted location in the Docker container.

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
[4] add_variant_position_alleles (optional) ──→ parsed positions/alleles
    ↓
[5] annotate_clinvar (optional) ──→ ClinVar clinical significance
    ↓
[6] annotate_gnomad (optional) ──→ gnomAD allele frequencies
    ↓
Output: fully annotated variants
```

### Step 1: Map Variants (Required)

**Purpose:** Normalize input variants (nucleotide or protein HGVS) to human-genome reference HGVS strings on GRCh38.

**Input columns:** `raw_hgvs_nt` (optional), `raw_hgvs_pro` (optional), `target_sequence` (required for sequence-based mapping)

**Output columns:** `mapped_hgvs_g`, `mapped_hgvs_c`, `mapped_hgvs_p`, `clingen_allele_id`, `mapped_errors`

**Command:**
```bash
src/scripts/run_map_variants.sh input.tsv output.tsv \
    --group-by gene_symbol \
    --raw-hgvs-nt raw_hgvs_nt \
    --raw-hgvs-pro raw_hgvs_pro
```

**Notes:**
- For transcript-referenced variants (e.g., `NM_000277.3:c.1218G>A`), the script normalizes through `dcd_mapping` to preserve accession handling
- For sequence-based variants without transcript accessions, same-gene rows are aligned together using BLAT + VRS
- Protein-only variants (no nucleotide HGVS) get protein mapping but c./g. fields will be empty
- See [BLAT Error 137 Retry Strategy](#blat-error-137-retry-strategy) section for handling memory issues

### Step 2: Reverse-Translate Protein Variants (Conditional)

**Purpose:** For protein-only variants mapped in step 1, reverse-translate them into multiple DNA (c./g.) candidates. This is essential when the same protein change can result from multiple DNA changes.

**Input columns:** `mapped_hgvs_p` (non-empty), `mapped_hgvs_c` and `mapped_hgvs_g` (must be empty/blank for this step to apply)

**Output columns:** Updates `mapped_hgvs_c` and `mapped_hgvs_g` with pipe-delimited candidates (e.g., `c.1218G>A|c.1221G>A`)

**Command:**
```bash
src/scripts/run_reverse_translate_protein_variants.sh output.tsv output_rt.tsv
```

**Notes:**
- Only processes rows where `mapped_hgvs_p` is present but `mapped_hgvs_c` and `mapped_hgvs_g` are empty
- Output candidates are pipe-delimited: `CA1||CA3` preserves position alignment (empty slots for unresolved candidates)
- DNA variants from step 1 pass through unchanged
- Queries transcript databases (UTA) to find all possible DNA backtranslations

### Step 3: Add DNA-Level ClinGen Allele IDs (Required for annotation)

**Purpose:** Resolve DNA-level ClinGen allele IDs for each reverse-translated candidate. For protein variants with multiple candidates, produces one ID per candidate. For DNA variants, reuses the existing ID from step 1.

**Input columns:** `mapped_hgvs_c`, `mapped_hgvs_g`, `clingen_allele_id` (from step 1)

**Output columns:** `dna_clingen_allele_id` (pipe-delimited to match candidates)

**Command:**
```bash
src/scripts/run_add_dna_clingen_allele_ids.sh output_rt.tsv output_clingen.tsv
```

**Notes:**
- For **DNA rows** (where `raw_hgvs_nt` is present): reuses existing `clingen_allele_id` if single candidate
- For **protein-only rows** (where `raw_hgvs_pro` is present): queries ClinGen for each candidate independently
- Falls back to g. (genomic) HGVS if c. (transcript) HGVS lookup fails
- Output is pipe-delimited to maintain alignment with reverse-translated candidates: `CA101||CA102` means 3 candidates, 2nd has no match
- This column becomes the primary key for downstream annotation steps

### Step 4: Add Parsed Position/Allele Columns (Optional)

**Purpose:** Parse HGVS strings into component fields (position, reference allele, alternate allele) for easier analysis and filtering.

**Input columns:** `mapped_hgvs_g`, `mapped_hgvs_c`, `mapped_hgvs_p`

**Output columns:** For each HGVS column, adds `<column>_start`, `<column>_stop`, `<column>_ref`, `<column>_alt`. Also adds `touches_intronic_region`, `spans_intron`

**Command:**
```bash
src/scripts/run_add_variant_position_alleles.sh output_clingen.tsv output_parsed.tsv
```

**Notes:**
- Useful for filtering by variant type or position
- Requires UTA database access for full ref-allele resolution (attempted but not strictly required)
- Safe to skip if you don't need parsed components

### Step 5: Annotate with ClinVar Data (Optional)

**Purpose:** Look up clinical significance, review status, and star ratings from ClinVar for DNA variants.

**Input columns:** `dna_clingen_allele_id` (from step 3)

**Output columns:** `clinvar.<YYYYMM>.clinical_significance`, `.review_status`, `.stars`, `.last_review_date`

**Command:**
```bash
src/scripts/run_annotate_clinvar.sh output_parsed.tsv output_clinvar.tsv \
    --clinvar-version 202601 \
    --cache-dir ./clinvar_cache
```

**Notes:**
- Uses pipe-delimited `dna_clingen_allele_id` candidates; tries each in order, returns first successful ClinVar match
- Downloads and caches the monthly ClinVar TSV from NCBI on first run
- Default namespace is `clinvar`; customize with `--namespace` flag
- Supports custom DNA ID column via `--dna-clingen-allele-id-col` if needed
- If all candidates fail to resolve, annotation columns are empty

### Step 6: Annotate with gnomAD Allele Frequencies (Optional)

**Purpose:** Look up allele frequencies and population-specific metrics from gnomAD using a local Hail table cache.

**Input columns:** `dna_clingen_allele_id` (from step 3)

**Output columns:** `gnomad.<VERSION>.minor_allele_frequency`, `.allele_frequency`, `.allele_count`, `.allele_number`, `.faf95_max`, `.faf95_max_ancestry`

**Command:**
```bash
src/scripts/run_annotate_gnomad.sh output_clinvar.tsv output_final.tsv \
    --gnomad-version v4.1 \
    --cache-dir ./gnomad_cache
```

**First-time setup (downloads and caches gnomAD Hail table):**
```bash
src/scripts/run_annotate_gnomad.sh output_clinvar.tsv output_final.tsv \
    --gnomad-version v4.1 \
    --cache-dir ./gnomad_cache \
    --download-only
```

**Notes:**
- Requires Java runtime and Hail library (installed via `gnomad` extra)
- First run downloads the source Hail table from GCP public data; this creates a local indexed cache
- Subsequent runs reuse the local cache (much faster)
- Default gnomAD version is `v4.1`; customize with `--gnomad-version` flag
- Custom table URI: `--gnomad-ht-uri gs://your-bucket/path/to/table.ht`
- Uses pipe-delimited candidates; tries each in order, returns first gnomAD match
- Supports custom DNA ID column via `--dna-clingen-allele-id-col` if needed
- Cache refresh: use `--refresh-cache` flag to re-download the source table

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
src/scripts/run_add_variant_position_alleles.sh variants_clingen.tsv variants_parsed.tsv

# Annotate with ClinVar (optional)
src/scripts/run_annotate_clinvar.sh variants_parsed.tsv variants_clinvar.tsv

# Annotate with gnomAD (optional, first time includes download)
src/scripts/run_annotate_gnomad.sh variants_clinvar.tsv variants_final.tsv
```

### Data Volume Management

ClinVar and gnomAD annotation steps use Docker volumes to persist caches:

```yaml
# In compose.yaml

volumes:
```
  variant-annotation-clinvar-cache:
  variant-annotation-gnomad-cache:

services:
  annotate-clinvar:
    volumes:
      - variant-annotation-clinvar-cache:/clinvar-cache

  annotate-gnomad:
    volumes:
      - variant-annotation-gnomad-cache:/gnomad-cache
```

This ensures caches survive container restarts and are shared across runs.

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
        │ add_variant_position   │   │
        │ _alleles               │   │
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
                       OUTPUT (annotated variants)
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

### Step 4: add_variant_position_alleles Issues

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

#### Step 4: add_variant_position_alleles

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
| `touches_intronic_region` | bool | True if transcript HGVS contains intron offsets |
| `spans_intron` | bool | True if transcript HGVS crosses an intron boundary |

#### Step 5: annotate_clinvar

| Column | Type | Description |
|--------|------|-------------|
| `clinvar.<YYYYMM>.clinical_significance` | string | ClinVar significance (e.g., `Pathogenic`, `Benign`, `Uncertain significance`) |
| `clinvar.<YYYYMM>.review_status` | string | Review status (e.g., `reviewed by expert panel`, `criteria provided`) |
| `clinvar.<YYYYMM>.stars` | int | Star rating: 0–4 (higher = more confident) |
| `clinvar.<YYYYMM>.last_review_date` | string | ISO 8601 date (e.g., `2023-06-15`) |

(Default namespace is `clinvar`; customize with `--namespace` flag)

#### Step 6: annotate_gnomad

| Column | Type | Description |
|--------|------|-------------|
| `gnomad.<VERSION>.minor_allele_frequency` | float | MAF = min(AF, 1-AF) |
| `gnomad.<VERSION>.allele_frequency` | float | Allele frequency in gnomAD cohort |
| `gnomad.<VERSION>.allele_count` | int | Count of alternate alleles in cohort |
| `gnomad.<VERSION>.allele_number` | int | Total alleles in cohort (2× sample count) |
| `gnomad.<VERSION>.faf95_max` | float | Filtering allele frequency (95% CI max) |
| `gnomad.<VERSION>.faf95_max_ancestry` | string | Ancestry group with highest FAF95 |

(Default version is `v4.1`; customize with `--gnomad-version` flag)

---

### Key Properties of Pipe-Delimited Columns

When a row has multiple DNA candidates (from reverse translation), the following columns are pipe-delimited:
- `mapped_hgvs_c`
- `mapped_hgvs_g`
- `dna_clingen_allele_id`

**Important:** Empty slots are preserved (e.g., `CA1||CA3` has 3 candidates, 2nd without a match).
This ensures positional alignment across all downstream annotations.

When annotating, the pipeline tries candidates in order and returns the first successful match.

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


