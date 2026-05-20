# clear_clingen_cache

`src/clear_clingen_cache.py` — delete ClinGen Allele Registry cache keys from Redis.

The pipeline caches ClinGen API responses in Redis under a namespace prefix (default: `clingen:v1`). Use this tool to invalidate the cache when you need to force fresh API lookups — for example after a ClinGen data update or when troubleshooting stale allele IDs.

Runs in Docker (uses the same `map-variants` container as `add_dna_clingen_allele_ids`).

---

## Usage

```bash
# Clear all keys under the default prefix (clingen:v1)
src/scripts/run_clear_clingen_cache.sh

# Clear keys under a specific prefix
src/scripts/run_clear_clingen_cache.sh --prefix clingen:v1
```

---

## CLI options

| Option | Default | Description |
|---|---|---|
| `--prefix PREFIX` | `CLINGEN_CACHE_PREFIX` env or `clingen:v1` | Redis key prefix to delete |
| `--log-level` | `INFO` | Logging verbosity |
