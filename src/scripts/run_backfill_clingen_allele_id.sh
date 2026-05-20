#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  cat <<'EOF'
Usage: src/scripts/run_backfill_clingen_allele_id.sh <input-file> <output-file> [options...]

Re-query the ClinGen Allele Registry for rows where clingen_allele_id is blank.
This fills the protein-level clingen_allele_id column produced by Step 1
(map_variants). It is NOT a substitute for Step 3 (add_dna_clingen_allele_ids),
which resolves DNA-level allele IDs for reverse-translated candidates.

Use this script when map_variants left some clingen_allele_id cells blank due
to transient API errors or rate-limiting during the initial run.

HGVS priority per row: mapped_hgvs_c > mapped_hgvs_g > mapped_hgvs_p
Apply only to Step-1 output (before reverse translation adds pipe-delimited values).

Examples:
  src/scripts/run_backfill_clingen_allele_id.sh input.tsv output.tsv
  src/scripts/run_backfill_clingen_allele_id.sh input.tsv output.tsv --concurrency 10
  src/scripts/run_backfill_clingen_allele_id.sh input.tsv output.tsv --skip 5000 --limit 1000

Key options (pass after output file; all are forwarded to the Python script):
  --clingen-allele-id COLNAME   Column to populate (default: clingen_allele_id)
  --concurrency N               Max concurrent ClinGen requests (default: 5)
  --max-retries N               Retries per request on failure (default: 3)
  --skip N                      Skip first N data rows
  --limit N                     Process at most N rows after --skip

Notes:
  - Paths are interpreted relative to /work in the container.
  - By default /work maps to ./data on the host.
  - Override mount directory with VARIANT_DATA_DIR=/absolute/path.
  - Images are reused by default for fast runs.
  - Add --rebuild-image to force rebuilding the tools image.
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

mapped_args=()
while [[ $# -gt 0 ]]; do
  mapped_args+=("$1")
  shift
done

cmd=(docker compose --profile tools run)
[[ -n "$compose_build_flag" ]] && cmd+=("$compose_build_flag")
[[ -n "$compose_no_cache_flag" ]] && cmd+=("$compose_no_cache_flag")
cmd+=(--rm --no-deps --entrypoint python map-variants
  -m src.backfill_clingen_allele_id "$input_in_container" "$output_in_container")
if [[ ${#mapped_args[@]} -gt 0 ]]; then
  cmd+=("${mapped_args[@]}")
fi
exec "${cmd[@]}"
