#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  cat <<'EOF'
Usage: src/scripts/run_annotate_gnomad.sh <input-file> <output-file> [annotate_gnomad options...]

Examples:
  src/scripts/run_annotate_gnomad.sh mapped.tsv annotated.tsv
  src/scripts/run_annotate_gnomad.sh mapped.tsv annotated.tsv --gnomad-version v4.1
  src/scripts/run_annotate_gnomad.sh mapped.tsv annotated.tsv --gnomad-ht-uri gs://.../gnomad.joint.v4.1.sites.ht

Notes:
  - Paths are interpreted relative to /work in the container.
  - By default /work maps to ./data on the host.
  - Override mount directory with VARIANT_DATA_DIR=/absolute/path.
  - By default annotate_gnomad reads DNA IDs from the dna_clingen_allele_id column.
  - Use --dna-clingen-allele-id-col to use a different input column name.
  - gnomAD Hail tables are cached in the variant-annotation-gnomad-cache Docker volume
    (mounted at /gnomad-cache). The GNOMAD_CACHE_DIR env var is set automatically.
  - Images are reused by default for fast runs.
  - Add --rebuild-image to force rebuilding the annotate-gnomad image.
  - Add --no-build-cache with --rebuild-image for a clean rebuild.
EOF
  exit 1
fi

input_path="$1"
output_path="$2"
shift 2

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
    --)
      shift
      break
      ;;
    *)
      break
      ;;
  esac
done

map_to_container_path() {
  local path="$1"
  local prefer_app_path="${2:-auto}"

  if [[ "$path" == /work/* || "$path" == /usr/src/app/* ]]; then
    printf '%s\n' "$path"
    return
  fi

  if [[ "$prefer_app_path" == "app" ]]; then
    printf '/usr/src/app/%s\n' "$path"
    return
  fi

  if [[ "$prefer_app_path" == "work" ]]; then
    printf '/work/%s\n' "$path"
    return
  fi

  if [[ -f "$path" ]]; then
    printf '/usr/src/app/%s\n' "$path"
  else
    printf '/work/%s\n' "$path"
  fi
}

input_in_container="$(map_to_container_path "$input_path")"

if [[ "$output_path" == /work/* || "$output_path" == /usr/src/app/* ]]; then
  output_in_container="$output_path"
elif [[ "$input_in_container" == /usr/src/app/* ]]; then
  output_in_container="/usr/src/app/$output_path"
else
  output_in_container="/work/$output_path"
fi

exec docker compose --profile tools run \
  ${compose_build_flag:+$compose_build_flag} ${compose_no_cache_flag:+$compose_no_cache_flag} --rm \
  annotate-gnomad "$input_in_container" "$output_in_container" "$@"
