#!/usr/bin/env bash
set -euo pipefail

cat_usage() {
  cat <<'EOF'
Usage: src/scripts/run_compare_cvfg_datasets.sh [options] [file_a [file_b [file_c]]]

All arguments are passed through to src.compare_cvfg_datasets.  File path
arguments are automatically mapped to /work (data directory) or /usr/src/app
(repository root) container paths as appropriate.

Positional arguments (all optional — Python defaults are used when omitted):
  file_a        Integrated variant effect dataset TSV
                  Default: data/cvfg/integrated_variant_effect_dataset.tsv
  file_b        CVFG variants flat TSV
                  Default: data/cvfg/v5/cvfg_variants.3.flat.tsv
  file_c        Dataset-to-score-set URN mapping TSV
                  Default: data/cvfg/dataset_to_score_set.tsv

Options (script-level, consumed before forwarding to Python):
  --output-dir <dir>     Output directory (forwarded to Python).  Path is
                         mapped to /work/<dir> unless already absolute.
  --rebuild-image        Force rebuild of the tools Docker image.
  --no-build-cache       (Use with --rebuild-image) disable layer cache.

All other flags (e.g. --no-check-ref-alleles) are forwarded to Python as-is.

Notes:
  - Paths are interpreted relative to /work in the container.
  - By default /work maps to ./data on the host.
  - Override mount directory with VARIANT_DATA_DIR=/absolute/path.
  - Images are reused by default for fast runs.
  - Requires UTA_DB_URL in settings/.env for the genomic ref-allele check.
EOF
}

if [[ $# -gt 0 && ( "$1" == "-h" || "$1" == "--help" ) ]]; then
  cat_usage
  exit 0
fi

compose_build_flag=""
compose_no_cache_flag=""

map_to_container_path() {
  local path="$1"
  local prefer_work="${2:-auto}"

  if [[ "$path" == /work/* || "$path" == /usr/src/app/* ]]; then
    printf '%s\n' "$path"
    return
  fi

  if [[ "$prefer_work" == "work" ]]; then
    printf '/work/%s\n' "$path"
    return
  fi

  if [[ "$prefer_work" == "app" ]]; then
    printf '/usr/src/app/%s\n' "$path"
    return
  fi

  # auto: existing files in the repo → /usr/src/app, everything else → /work
  if [[ -f "$path" ]]; then
    printf '/usr/src/app/%s\n' "$path"
  else
    printf '/work/%s\n' "$path"
  fi
}

# ── Parse script-level flags and collect Python args ──────────────────────────
python_args=()
positional_count=0

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
    --output-dir)
      [[ $# -lt 2 ]] && { echo "Error: --output-dir requires a value"; exit 1; }
      python_args+=("--output-dir" "$(map_to_container_path "$2" work)")
      shift 2
      ;;
    --output-dir=*)
      python_args+=("--output-dir" "$(map_to_container_path "${1#*=}" work)")
      shift
      ;;
    --)
      shift
      # Everything after -- is forwarded verbatim
      while [[ $# -gt 0 ]]; do
        python_args+=("$1")
        shift
      done
      ;;
    -*)
      # Forward unknown flags verbatim
      python_args+=("$1")
      shift
      ;;
    *)
      # Positional: file_a, file_b, file_c (in order)
      positional_count=$((positional_count + 1))
      if [[ $positional_count -le 3 ]]; then
        python_args+=("$(map_to_container_path "$1")")
      else
        python_args+=("$1")
      fi
      shift
      ;;
  esac
done

# ── Build and run the Docker Compose command ───────────────────────────────────
cmd=(docker compose --profile tools run)
[[ -n "$compose_build_flag" ]] && cmd+=("$compose_build_flag")
[[ -n "$compose_no_cache_flag" ]] && cmd+=("$compose_no_cache_flag")
cmd+=(--rm compare-cvfg-datasets)
if [[ ${#python_args[@]} -gt 0 ]]; then
  cmd+=("${python_args[@]}")
fi
exec "${cmd[@]}"
