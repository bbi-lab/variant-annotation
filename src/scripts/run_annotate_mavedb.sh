#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  cat <<'EOF'
Usage: src/scripts/run_annotate_mavedb.sh <input-file> <output-file> [annotate_mavedb options...]

Examples:
  src/scripts/run_annotate_mavedb.sh variants.tsv annotated.tsv
  src/scripts/run_annotate_mavedb.sh variants.tsv annotated.tsv \
    --mavedb-api-url https://api.mavedb.org
  src/scripts/run_annotate_mavedb.sh /work/variants.tsv /work/annotated.tsv \
    --variant-urn-col variant_urn \
    --score-col score

Notes:
  - Input/output paths are interpreted relative to /work in the container.
  - By default /work maps to ./data on the host.
  - Override mount directory with VARIANT_DATA_DIR=/absolute/path.
  - The MaveDB API URL defaults to https://api.mavedb.org.  Set MAVEDB_API_URL
    in settings/.env to override (e.g. to point at a local instance).
  - For each unique score set URN found in the variant_urn column the script
    calls the MaveDB REST API once to retrieve calibrations.  For class-based
    calibrations an additional per-calibration API call is made to fetch variant
    class assignments.
  - Images are reused by default for fast runs.
  - Add --rebuild-image to force rebuilding the image.
  - Add --no-build-cache with --rebuild-image for a clean rebuild.

Output columns added:
  mavedb.primary_calibration.urn
  mavedb.primary_calibration.name
  mavedb.primary_calibration.functional_class
  mavedb.investigator_provided_calibration.urn
  mavedb.investigator_provided_calibration.name
  mavedb.investigator_provided_calibration.functional_class
EOF
  exit 1
fi

input_path="$1"
output_path="$2"
shift 2

compose_build_flag=""
compose_no_cache_flag=""

# Wrapper-only options (must appear before annotate_mavedb options).
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
  annotate-mavedb "$input_in_container" "$output_in_container" "$@"
