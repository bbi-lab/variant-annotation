# annotate_gnomad

`src/annotate_gnomad.py` — Step 6 of the variant-annotation pipeline (optional).

Annotates each variant row with gnomAD allele frequency metrics. Supports two execution backends (`--execution-mode hail` or `athena`) and two variant lookup strategies (`--lookup-mode coordinates` or `caid`).

---

## Output columns

Output column names follow `<namespace>.<version>.<field>` (namespace defaults to `gnomad`, version defaults to `v4.1`).

### Standard columns (always emitted)

| Column | Description |
|---|---|
| `gnomad.<V>.minor_allele_frequency` | Minor allele frequency (= `min(AF, 1-AF)`) |
| `gnomad.<V>.allele_frequency` | Allele frequency (`AC / AN`) |
| `gnomad.<V>.allele_count` | Allele count (`AC`) |
| `gnomad.<V>.allele_number` | Allele number (`AN`) |
| `gnomad.<V>.faf95_max` | Maximum filtering allele frequency across ancestry groups (95% CI) |
| `gnomad.<V>.faf95_max_ancestry` | Ancestry group that produced `faf95_max` |
| `gnomad.<V>.filters` | Combined gnomAD QC filters (pipe-delimited when multiple; empty = PASS) |
| `gnomad.<V>.exome_filters` | Exome-callset QC filters (empty = passed or not called in exomes) |
| `gnomad.<V>.genome_filters` | Genome-callset QC filters (empty = passed or not called in genomes) |
| `gnomad.<V>.gene_symbols` | Gene symbols from `vep.worst_csq_by_gene_canonical` (pipe-delimited) |

All standard columns are pipe-delimited and candidate-aligned when `dna_clingen_allele_id` or the coordinate columns are themselves pipe-delimited (i.e. for rows with multiple reverse-translated DNA candidates from step 2).

### Optional histogram columns (Hail mode only)

When `--age-histograms` or `--allele-balance-histograms` is specified, one group of columns is added per requested histogram. The column names are derived at runtime from the bin edges found in the first matched gnomAD record, and therefore reflect the actual gnomAD bin boundaries rather than hard-coded values.

**Column naming convention:**

```
gnomad.<V>.<histogram_name>.n_smaller
gnomad.<V>.<histogram_name>.bin_1_<lo>_<hi>
gnomad.<V>.<histogram_name>.bin_2_<lo>_<hi>
...
gnomad.<V>.<histogram_name>.bin_N_<lo>_<hi>
gnomad.<V>.<histogram_name>.n_larger
```

Where `<lo>` and `<hi>` are the lower and upper bin edges formatted without trailing zeros (e.g. `17.5`, `25`, `0.05`).

**Age distribution histograms** (`--age-histograms`):

| `<histogram_name>` | Callset | Genotype |
|---|---|---|
| `age_hist_exome_het` | Exome | Heterozygous |
| `age_hist_exome_hom` | Exome | Homozygous |
| `age_hist_genome_het` | Genome | Heterozygous |
| `age_hist_genome_hom` | Genome | Homozygous |
| `age_hist_joint_het` | Joint | Heterozygous |
| `age_hist_joint_hom` | Joint | Homozygous |

For gnomAD v4.1 the age bins span 17.5 – 80 years in irregular intervals, yielding 12 bins + `n_smaller` + `n_larger` = 14 values per histogram.

**Allele balance histograms** (`--allele-balance-histograms`):

| `<histogram_name>` | Callset | Filter level |
|---|---|---|
| `ab_hist_exome_adj` | Exome | Adjusted (genotype-filtered) |
| `ab_hist_exome_raw` | Exome | Raw (unfiltered) |
| `ab_hist_genome_adj` | Genome | Adjusted |
| `ab_hist_genome_raw` | Genome | Raw |
| `ab_hist_joint_adj` | Joint | Adjusted |
| `ab_hist_joint_raw` | Joint | Raw |

For gnomAD v4.1 the allele balance bins span 0 – 1 in 20 equal-width intervals, yielding 20 bins + `n_smaller` + `n_larger` = 22 values per histogram.

Histogram columns are also pipe-delimited and candidate-aligned, matching the standard columns.

---

## Execution modes

### Hail mode (default)

Downloads the gnomAD Hail table from a source URI (GCS or local filesystem), builds a compact local indexed cache, and queries that cache for each run. Requires Java and the `hail` Python library.

**First run** (building the cache — standard columns only):
```bash
src/scripts/run_annotate_gnomad.sh /dev/null /dev/null \
    --gnomad-version v4.1 \
    --download-only \
    --refresh-cache \
    --gnomad-ht-uri gs://gcp-public-data--gnomad/release/4.1/ht/joint/gnomad.joint.v4.1.sites.ht
```

The first run writes a compact local Hail table to `--cache-dir`. This is a one-time step and typically takes hours for a full gnomAD joint sites table (GCS read + Parquet write). The cache is a Hail keyed table (`gnomad_{version}_indexed.ht`) stored as a directory. Subsequent runs use the cache directly in seconds.

**First run with all histogram columns included in the cache:**
```bash
src/scripts/run_annotate_gnomad.sh /dev/null /dev/null \
    --gnomad-version v4.1 \
    --download-only \
    --refresh-cache \
    --gnomad-ht-uri gs://gcp-public-data--gnomad/release/4.1/ht/joint/gnomad.joint.v4.1.sites.ht \
    --age-histograms exome,genome,joint \
    --allele-balance-histograms exome,genome,joint
```

Histogram field names passed at cache-build time determine which histogram arrays are stored in the local Hail table. If you request histograms on an annotation run against a cache built without them, a warning is logged and those columns are omitted. Use `--refresh-cache` to rebuild the cache with the new fields.

To cache a local gnomAD HT copy instead of downloading from GCS:
```bash
src/scripts/run_annotate_gnomad.sh /dev/null /dev/null \
    --gnomad-version v4.1 \
    --download-only --refresh-cache \
    --gnomad-ht-uri /work/gnomAD/gnomad.joint.v4.1.sites.ht
```

**Annotation run (standard columns):**
```bash
src/scripts/run_annotate_gnomad.sh input.tsv output.tsv \
    --gnomad-version v4.1
```

**Annotation run with all histograms enabled:**
```bash
src/scripts/run_annotate_gnomad.sh input.tsv output.tsv \
    --gnomad-version v4.1 \
    --age-histograms exome,genome,joint \
    --allele-balance-histograms exome,genome,joint
```

**Selective histogram output (exome age distribution only):**
```bash
src/scripts/run_annotate_gnomad.sh input.tsv output.tsv \
    --gnomad-version v4.1 \
    --age-histograms exome
```

### Athena mode

Queries gnomAD via AWS Athena. No local cache or Hail installation is required. Intended for environments where gnomAD is already registered as an Athena external table (e.g. a shared MaveDB worker environment).

```bash
src/scripts/run_annotate_gnomad.sh input.tsv output.tsv \
    --execution-mode athena \
    --gnomad-version v4.1 \
    --athena-output-location s3://your-bucket/athena-results/
```

Input rows are processed in batches; results are written in input order after each batch.

**Important:** QC filtering (`--require-pass`, `--callset-pass-filter`) is **not available** in Athena mode. The Athena query only fetches `AC`, `AN`, `faf95_max`, and `faf95_max_ancestry` — it does not retrieve filter columns. The `filters`, `exome_filters`, and `genome_filters` output columns are always empty in Athena mode, and the filter flags are silently no-ops. Use Hail mode if you need QC filtering.

---

## Lookup modes

### `--lookup-mode coordinates` (default)

Builds a `"chrN:pos:ref:alt"` key from the VCF-style coordinate columns (`mapped_hgvs_g_chromosome`, `mapped_hgvs_g_start`, `mapped_hgvs_g_ref`, `mapped_hgvs_g_alt`) and looks up gnomAD records by that key.

This is the preferred mode because many gnomAD records have no CAID, so CAID-based lookup misses them entirely.

### `--lookup-mode caid`

Uses `dna_clingen_allele_id` values directly. When the local gnomAD cache is keyed by CAID (case 1 below), this is a direct index lookup. When the cache is coordinate-keyed, the CAID must first be resolved to coordinates.

---

## CAID-to-record resolution (Hail mode, `--lookup-mode caid`)

Three strategies are tried in order depending on the local cache structure:

| Case | Condition | Strategy |
|---|---|---|
| 1 | Source table has a `caid` field | Direct CAID index lookup in the local cache |
| 2 | Coordinate columns present in the input file | Pre-computed `CAID → chrN:pos:ref:alt` mapping from input columns |
| 3 | Fallback | ClinGen Allele Registry API → GRCh38 coordinates → cache lookup |

Case 3 makes live HTTP requests to the ClinGen API and can be slow for large batches.

### CAID normalisation

gnomAD stores CAIDs without leading zeroes (e.g. `CA25094`), while the ClinGen Allele Registry returns them zero-padded (e.g. `CA025094`). The script strips leading zeroes from the numeric suffix before any lookup.

---

## QC filtering (Hail mode only)

Two independent filtering flags are applied at annotation time:

### `--require-pass`

Excludes variants whose gnomAD **combined** `filters` field is non-empty. The combined field aggregates QC filters across all callsets (exomes + genomes). An empty `filters` field means the variant passed all gnomAD QC — equivalent to `FILTER == PASS` in a joint VCF.

Common gnomAD filter values: `AC0` (allele count zero), `AS_VQSR` (VQSR failure), `InbreedingCoeff`.

### `--callset-pass-filter {none|any|all}`

Filters on per-callset QC independently:

| Value | Behaviour |
|---|---|
| `none` (default) | No callset filtering; all matched records are annotated |
| `any` | Only annotate variants that passed QC in **at least one** callset (exome OR genome `filters` field is empty) |
| `all` | Only annotate variants that passed QC in **both** callsets (exome AND genome `filters` fields are both empty) |

An empty `exome_filters` or `genome_filters` field means the variant either passed in that callset or was not called in it. Non-empty means one or more filters were applied (e.g. `AC0|AS_VQSR`).

`--callset-pass-filter any/all` requires that the gnomAD cache was built from a table with separate `exome.filters` / `genome.filters` fields (present in the gnomAD v4.1 joint sites table). If the cache was built from a table without those fields, an error is raised and you should use `--require-pass` instead.

`--require-pass` and `--callset-pass-filter` can be combined; a variant must satisfy **both** conditions.

A variant with no gnomAD record is always left unannotated regardless of filter settings.

---

## Cache management

| Flag | Effect |
|---|---|
| `--download-only` | Build or refresh the local cache without annotating any rows. Ignored in Athena mode. |
| `--refresh-cache` | Force rebuilding the local cache even if one already exists. Stale Hail shuffle data is cleaned automatically. Ignored in Athena mode. |
| `--cache-progress-every-seconds N` | Log cache preparation progress every N seconds (default 300; 0 disables). Useful for the multi-hour initial cache build. |

The local cache directory defaults to `$GNOMAD_CACHE_DIR` or `/tmp/gnomad_cache`. In Docker Compose the cache is mapped to a named volume. The default source URI is the gnomAD v4.1 joint sites GCS bucket, configurable via `GNOMAD_HT_URI` or `--gnomad-ht-uri`.

### Gene-level cache filtering

Pass `--genes BRCA1,BRCA2` when building the cache to retain only rows whose `vep.worst_csq_by_gene_canonical` contains at least one matching gene symbol. Reduces cache size by ~100× for targeted analyses. Requires `--refresh-cache` (or a cache miss) to take effect; an existing cache is not re-filtered.

### Histogram fields in the cache

Histogram arrays (`--age-histograms`, `--allele-balance-histograms`) are stored as array fields in the local cache. They are only present when explicitly requested at cache-build time. If you run an annotation with histogram flags against an older cache, the tool will log a warning per missing histogram name and omit those columns from output. Rebuild the cache with `--refresh-cache` (and the histogram flags) to include them.

The cache is a superset: you can request a subset of histograms on any given annotation run, as long as all requested fields were included when the cache was built.

---

## CLI options

### Common options

| Option | Default | Description |
|---|---|---|
| `--execution-mode {hail,athena}` | `hail` | Execution backend |
| `--lookup-mode {coordinates,caid}` | `coordinates` | Variant lookup strategy |
| `--gnomad-version VERSION` | `v4.1` | Version label for output column names |
| `--gnomad-namespace NS` | `gnomad` | Prefix for output column names |
| `--dna-clingen-allele-id-col COL` | `dna_clingen_allele_id` | CAID input column (`--lookup-mode caid`) |
| `--skip N` | `0` | Skip the first N data rows |
| `--limit N` | no limit | Stop after N rows |
| `--log-level` | `INFO` | Logging verbosity |
| `--csv-field-size-limit BYTES` | system default | Increase for large HGVS fields |

### QC filtering (Hail mode only)

| Option | Default | Description |
|---|---|---|
| `--require-pass` | off | Exclude variants with non-empty combined gnomAD `filters` field |
| `--callset-pass-filter {none,any,all}` | `none` | Per-callset QC filter level |

### Coordinate lookup columns (`--lookup-mode coordinates`)

| Option | Default |
|---|---|
| `--coord-chromosome-col` | `mapped_hgvs_g_chromosome` |
| `--coord-pos-col` | `mapped_hgvs_g_start` |
| `--coord-ref-col` | `mapped_hgvs_g_ref` |
| `--coord-alt-col` | `mapped_hgvs_g_alt` |

### Hail-mode options

| Option | Default | Description |
|---|---|---|
| `--gnomad-ht-uri URI` | `$GNOMAD_HT_URI` or GCS public path | Source gnomAD Hail table |
| `--cache-dir DIR` | `$GNOMAD_CACHE_DIR` or `/tmp/gnomad_cache` | Local cache directory |
| `--download-only` | off | Build cache without annotating |
| `--refresh-cache` | off | Rebuild cache even if present |
| `--genes GENE,...` | none | Gene-level cache filter |
| `--cache-progress-every-seconds N` | `300` | Cache build heartbeat interval |
| `--age-histograms CALLSET,...` | none | Include age distribution histogram columns; comma-separated from `exome`, `genome`, `joint` |
| `--allele-balance-histograms CALLSET,...` | none | Include allele balance histogram columns; comma-separated from `exome`, `genome`, `joint` |

### Athena-mode options

| Option | Env variable | Default | Description |
|---|---|---|---|
| `--athena-database` | `GNOMAD_ATHENA_DATABASE` | `gnomad` | Athena database name |
| `--athena-table` | `GNOMAD_ATHENA_TABLE` | derived from `--gnomad-version` | Athena table name |
| `--athena-output-location` | `GNOMAD_ATHENA_OUTPUT_LOCATION` | — | S3 URI for Athena query results (**required**) |
| `--athena-workgroup` | `GNOMAD_ATHENA_WORKGROUP` | — | Optional Athena workgroup |
| `--athena-region` | `GNOMAD_ATHENA_REGION` / `AWS_REGION` | — | AWS region for Athena client |
| `--athena-max-caids-per-query N` | — | `16250` | Max CAIDs per Athena `IN` clause |
| `--athena-row-batch-size N` | `GNOMAD_ATHENA_ROW_BATCH_SIZE` | `1000` | Input rows per batch |
| `--athena-poll-seconds N` | — | `5` | Athena query polling interval |
| `--athena-coord-query-strategy` | — | `by-chromosome` | SQL strategy for coordinate lookups: `by-chromosome` (partition-pruned per-chromosome queries) or `or-of-ands` (exact four-column match per batch) |

---

## Dependencies

### Hail mode

- **Java runtime** (JVM) — required by Apache Spark / Hail.
- **`hail` Python library** — install via `pip install -e '.[gnomad]'` or the Docker image.
- **GCS connector** (auto-configured) — required only for `gs://` source URIs. The Docker image configures `spark.hadoop.fs.gs.*` and downloads the connector JAR automatically. Set `HAIL_GCS_AUTH_TYPE=UNAUTHENTICATED` for public gnomAD buckets.
- **Spark local directory** — Spark shuffle data defaults to `$SPARK_LOCAL_DIR` (or the Hail tmp dir). Point it at the Docker volume (not the overlay filesystem) to avoid disk-full hangs during the cache write.

### Athena mode

- **`boto3`** — included in project dependencies. AWS credentials must be available via the standard boto3 credential chain (environment variables, IAM instance profile, `~/.aws/credentials`, etc.).

---

## Troubleshooting

**Cache build hangs or runs out of disk space**

Set `SPARK_LOCAL_DIR` to a path inside the Docker volume mount. Spark shuffle data during the cache write stage can reach 26+ GiB. The Docker image maps `$GNOMAD_CACHE_DIR` to a named volume; setting `SPARK_LOCAL_DIR=$GNOMAD_CACHE_DIR/spark-tmp` keeps shuffle data off the container overlay filesystem.

**`UnsupportedFileSystemException: scheme "gs"`**

The GCS connector was not found. Rebuild the Docker image (`--rebuild-image`). If the problem persists in your environment, download the gnomAD HT locally and use a local `--gnomad-ht-uri`.

**`No valid credential configuration discovered`**

Set `HAIL_GCS_AUTH_TYPE=UNAUTHENTICATED` for the public gnomAD GCS buckets. This is the default in the Docker image.

**`--callset-pass-filter any/all` raises an error**

The gnomAD cache was built from a source table without separate `exome.filters` / `genome.filters` fields (e.g. the browser sites HT). Use `--require-pass` instead, or rebuild the cache from the gnomAD joint sites HT.

**Histogram columns missing from output despite passing `--age-histograms` / `--allele-balance-histograms`**

The local cache was built without those histogram fields. The log will contain a warning such as:
```
Histogram field(s) age_hist_exome_het, age_hist_exome_hom not found in gnomAD cache; use --refresh-cache to rebuild with histogram support
```
Rebuild the cache with `--refresh-cache` and the same histogram flags, then re-run annotation.

**Inconsistent bin edges error**

```
ValueError: Inconsistent bin edges for histogram 'age_hist_exome_het' at key ...
```

The gnomAD source table contains variant rows with different bin-edge arrays for the same histogram field. This should not occur in standard gnomAD release tables but may happen with custom or merged tables. If you encounter it, open an issue with the source table URI.

**Many empty frequency columns despite valid `dna_clingen_allele_id`**

- Check `--lookup-mode`. The default `coordinates` mode uses VCF coordinate columns; if those are blank, no match is possible. Confirm `mapped_hgvs_g_chromosome`, `mapped_hgvs_g_start`, `mapped_hgvs_g_ref`, `mapped_hgvs_g_alt` are populated from step 4.
- With `--lookup-mode caid`, check for CAID normalisation (leading-zero stripping) and verify the gnomAD cache was keyed by CAID.
- Intronic variants and some indels are not covered by gnomAD; empty results are expected for those.
