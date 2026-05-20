# gnomad_hail_shell

`src/scripts/run_gnomad_hail_shell.sh` — launch an interactive bash shell in the Hail-enabled gnomAD container.

Opens a `bash` session inside the `annotate-gnomad` Docker container with the same mounts used during pipeline annotation runs. Useful for ad-hoc inspection of gnomAD Hail tables, testing queries, or diagnosing annotation issues.

---

## Mounts

| Container path | Host path |
|---|---|
| `/work` | `./data` (or `$VARIANT_DATA_DIR`) |
| `/gnomad-cache` | `variant-annotation-gnomad-cache` Docker volume |

Hail environment variables (`HAIL_SPARK_JARS`, `HAIL_GCS_CONNECTOR_JAR`) are pre-configured in the image.

---

## Usage

```bash
src/scripts/run_gnomad_hail_shell.sh

# Rebuild the image first
src/scripts/run_gnomad_hail_shell.sh --rebuild-image
```

| Option | Description |
|---|---|
| `--rebuild-image` | Rebuild the `annotate-gnomad-hail` image before starting the shell |
| `--no-build-cache` | Disable Docker layer cache (use with `--rebuild-image`) |
