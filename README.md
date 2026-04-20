# variant-annotation

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
