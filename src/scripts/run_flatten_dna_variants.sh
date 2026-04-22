#!/bin/bash
#
# Wrapper script for flattening pipe-delimited DNA variants to one row per DNA candidate.
#
# Converts annotated variant TSV files (with pipe-delimited DNA candidates) to a flattened
# format where each DNA candidate has its own row. Protein-only variants without reverse
# translations are dropped.
#
# Usage:
#   src/scripts/run_flatten_dna_variants.sh <input.tsv> <output.tsv> [OPTIONS]
#
# Options:
#   --rebuild-image              Rebuild Docker image before running
#   --no-build-cache             Build without Docker cache
#   --dna-variant-columns COLS   Comma-separated list of columns to expand (auto-detect if omitted)
#
# Examples:
#   # Flatten with auto-detected columns
#   src/scripts/run_flatten_dna_variants.sh annotated.tsv dna_variants.tsv
#
#   # Explicitly specify columns to expand
#   src/scripts/run_flatten_dna_variants.sh annotated.tsv dna_variants.tsv \
#     --dna-variant-columns mapped_hgvs_g,mapped_hgvs_c,dna_clingen_allele_id
#
#   # Rebuild image before running
#   src/scripts/run_flatten_dna_variants.sh annotated.tsv dna_variants.tsv --rebuild-image
#

set -euo pipefail

if [[ $# -lt 2 ]]; then
  cat <<'EOF'
Usage: src/scripts/run_flatten_dna_variants.sh <input-file> <output-file> [flatten_dna_variants options...]

Examples:
  src/scripts/run_flatten_dna_variants.sh input.tsv output.tsv
  src/scripts/run_flatten_dna_variants.sh input.tsv output.tsv --dna-variant-columns col1,col2

Notes:
  - Paths are interpreted relative to /work in the container.
  - By default /work maps to ./data on the host.
  - Override mount directory with VARIANT_DATA_DIR=/absolute/path.
  - Images are reused by default for fast runs.
  - Add --rebuild-image to force rebuilding the flatten-dna-variants image.
  - Add --no-build-cache with --rebuild-image for a clean rebuild.
EOF
  exit 1
fi

input_path="$1"
output_path="$2"
shift 2

compose_build_flag=""
compose_no_cache_flag=""

# Wrapper-only options (must appear before flatten_dna_variants options).
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

mapped_args=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --dna-variant-columns)
      mapped_args+=("$1")
      if [[ $# -lt 2 ]]; then
        echo "Error: missing value for $1" >&2
        exit 1
      fi
      mapped_args+=("$2")
      shift 2
      ;;
    *)
      mapped_args+=("$1")
      shift
      ;;
  esac
done

exec docker compose --profile tools run \
  ${compose_build_flag:+$compose_build_flag} ${compose_no_cache_flag:+$compose_no_cache_flag} --rm \
  flatten-dna-variants "$input_in_container" "$output_in_container" "${mapped_args[@]}"
