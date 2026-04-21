# variant-annotation

## License

This project is licensed under the GNU Affero General Public License v3.0 or
later. See the [LICENSE](LICENSE) file.

The AGPL choice reflects that parts of the variant-mapping pipeline are
intentionally modeled on established MaveDB workflow behavior and may include
adapted implementation ideas or logic in that area of the codebase.

## Prerequisites

- [Docker](https://docs.docker.com/get-docker/) and [Docker Compose](https://docs.docker.com/compose/install/)
- Python 3.11+

## Development setup

The development environment runs entirely in Docker. All services (app, cdot, database, seqrepo, UTA) are defined in `docker-compose-dev.yml`.

For TSV processing, use the dedicated `map-variants` service. It mounts a host folder at `/work` in the container:

- default host folder: `./data`
- optional override: set `VARIANT_DATA_DIR=/absolute/path/to/your/files`
- on Apple Silicon, `map-variants` runs as `linux/amd64` so UCSC BLAT works
- the Docker image installs `polars[rtcompat]` to avoid AVX-related crashes under emulation
- class-1 (accession-referenced) nucleotide HGVS variants are normalized through `dcd_mapping` (matching MaveDB behavior) before ClinGen extraction, which preserves accession/version handling for ENST inputs

### 1. Configure environment

Copy the example env file and fill in the required values:

```bash
cp settings/.env.example settings/.env
cp settings/.env.example settings/.env.dev
```

### 2. Build and start services

Before the first startup, fetch the UTA database dump once:

```bash
src/scripts/fetch_uta_dump.sh
```

Then start services:

```bash
docker compose -f docker-compose-dev.yml up --build
```

The `cdot` service is built directly from source at [github.com/VariantEffect/cdot_rest](https://github.com/VariantEffect/cdot_rest) — no manual clone required.

### 3. Open a shell in the app container

```bash
docker compose -f docker-compose-dev.yml exec app bash
```

### 4. Run map_variants on a TSV from the host

Put your input file in `./data` (or your `VARIANT_DATA_DIR` folder), then run:

```bash
docker compose -f docker-compose-dev.yml --profile tools run --rm map-variants /work/input.tsv /work/output.tsv
```

Shortcut wrapper script:

```bash
src/scripts/run_map_variants.sh input.tsv output.tsv
```

If `input.tsv` exists in the repository working tree, the wrapper automatically
uses the project bind mount (`/usr/src/app/...`). Otherwise it defaults to
the `/work` mount backed by `./data` (or `VARIANT_DATA_DIR`).

Example with options:

```bash
docker compose -f docker-compose-dev.yml --profile tools run --rm map-variants \
	/work/input.tsv \
	/work/output.tsv \
	--group-by gene_symbol \
	--raw-hgvs-nt raw_hgvs_nt \
	--raw-hgvs-pro raw_hgvs_pro
```

Equivalent with the wrapper script:

```bash
src/scripts/run_map_variants.sh input.tsv output.tsv \
	--group-by gene_symbol \
	--raw-hgvs-nt raw_hgvs_nt \
	--raw-hgvs-pro raw_hgvs_pro
```

Reuse existing mapped results to avoid reprocessing matching variants:

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

## Troubleshooting

### BLAT Error 137 (SIGKILL / Out of Memory)

**Symptom:** Mapping fails with error message: "BLAT process returned error code 137"

**Root Cause:** The BLAT subprocess was killed by the operating system due to out-of-memory
(OOM) conditions. This typically occurs when processing very large groups of variants in a
single BLAT run, as each variant pair requires significant memory for alignment.

**Solution:** The `map-variants` command automatically enables error-driven retry with
progressive chunking by default. No manual intervention is usually needed. In cases where
error 137 still persists:

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
   Edit `docker-compose-dev.yml` and add memory limits under the `map-variants` service:
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
