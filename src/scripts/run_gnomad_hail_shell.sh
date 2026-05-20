#!/usr/bin/env bash
set -euo pipefail

# Start an interactive bash shell in the Hail-enabled gnomAD container.
# The gnomad cache volume and data directory are mounted as usual.
#
# Usage:
#   src/scripts/run_gnomad_hail_shell.sh [--rebuild-image] [--no-build-cache]
#
# Notes:
#   - /work maps to ./data (or $VARIANT_DATA_DIR on the host)
#   - /gnomad-cache maps to the variant-annotation-gnomad-cache Docker volume
#   - HAIL_SPARK_JARS and HAIL_GCS_CONNECTOR_JAR are pre-configured in the image
#   - Add --rebuild-image to rebuild the annotate-gnomad-hail image first

compose_build_flag=""
compose_no_cache_flag=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --rebuild-image)
      compose_build_flag="--build"
      shift
      ;;
    --no-build-cache)
      compose_no_cache_flag="--no-cache"
      shift
      ;;
    *)
      echo "Unknown option: $1" >&2
      exit 1
      ;;
  esac
done

exec docker compose --profile tools run \
  ${compose_build_flag:+$compose_build_flag} ${compose_no_cache_flag:+$compose_no_cache_flag} \
  --rm -it --entrypoint bash annotate-gnomad
