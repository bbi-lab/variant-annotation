#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  cat <<'EOF'
Usage: src/scripts/run_map_variants.sh <input-file> <output-file> [map_variants options...]

Examples:
  src/scripts/run_map_variants.sh input.tsv output.tsv
  src/scripts/run_map_variants.sh input.tsv output.tsv --group-by gene_symbol
  src/scripts/run_map_variants.sh input.tsv output.tsv --targets-file data/targets.tsv

Notes:
  - Paths are interpreted relative to /work in the container.
  - By default /work maps to ./data on the host.
  - Override mount directory with VARIANT_DATA_DIR=/absolute/path.
  - Images are reused by default for fast runs.
  - Add --rebuild-image to force rebuilding the map-variants image.
  - Add --no-build-cache with --rebuild-image for a clean rebuild.
  - --targets-file and --merge-existing paths are mapped to the container automatically.
EOF
  exit 1
fi

input_path="$1"
output_path="$2"
shift 2

compose_build_flag=""
compose_no_cache_flag=""

# Wrapper-only options (must appear before map_variants options).
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
    --merge-existing)
      if [[ $# -lt 2 ]]; then
        echo "error: --merge-existing requires a file argument" >&2
        exit 2
      fi
      mapped_args+=("--merge-existing")
      if [[ "$input_in_container" == /usr/src/app/* ]]; then
        mapped_args+=("$(map_to_container_path "$2" app)")
      else
        mapped_args+=("$(map_to_container_path "$2")")
      fi
      shift 2
      ;;
    --merge-existing=*)
      merge_path="${1#*=}"
      if [[ "$input_in_container" == /usr/src/app/* ]]; then
        mapped_args+=("--merge-existing=$(map_to_container_path "$merge_path" app)")
      else
        mapped_args+=("--merge-existing=$(map_to_container_path "$merge_path")")
      fi
      shift
      ;;
    --targets-file)
      if [[ $# -lt 2 ]]; then
        echo "error: --targets-file requires a file argument" >&2
        exit 2
      fi
      mapped_args+=("--targets-file")
      if [[ "$input_in_container" == /usr/src/app/* ]]; then
        mapped_args+=("$(map_to_container_path "$2" app)")
      else
        mapped_args+=("$(map_to_container_path "$2")")
      fi
      shift 2
      ;;
    --targets-file=*)
      targets_path="${1#*=}"
      if [[ "$input_in_container" == /usr/src/app/* ]]; then
        mapped_args+=("--targets-file=$(map_to_container_path "$targets_path" app)")
      else
        mapped_args+=("--targets-file=$(map_to_container_path "$targets_path")")
      fi
      shift
      ;;
    *)
      mapped_args+=("$1")
      shift
      ;;
  esac
done

exec docker compose --profile tools run \
  ${compose_build_flag:+$compose_build_flag} ${compose_no_cache_flag:+$compose_no_cache_flag} --rm \
  map-variants "$input_in_container" "$output_in_container" "${mapped_args[@]}"
