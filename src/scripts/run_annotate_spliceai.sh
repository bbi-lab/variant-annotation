#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  cat <<'EOF'
Usage: src/scripts/run_annotate_spliceai.sh <input-file> <output-file> [annotate_spliceai options...]

Examples:
  src/scripts/run_annotate_spliceai.sh mapped.tsv annotated.tsv
  src/scripts/run_annotate_spliceai.sh mapped.tsv annotated.tsv --mode precomputed
  src/scripts/run_annotate_spliceai.sh mapped.tsv annotated.tsv \
    --precomputed-snv-vcf spliceai_scores.masked.snv.hg38.vcf.gz \
    --precomputed-indel-vcf spliceai_scores.masked.indel.hg38.vcf.gz

Notes:
  - Paths are interpreted relative to /work in the container.
  - By default /work maps to ./data on the host.
  - Override mount directory with VARIANT_DATA_DIR=/absolute/path.
  - Precomputed SpliceAI files are copied/indexed into the variant-annotation-spliceai-cache
    Docker volume (mounted at /spliceai-cache).
  - If .tbi indexes are missing, they are generated automatically.
  - Images are reused by default for fast runs.
  - Add --rebuild-image to force rebuilding the annotate-spliceai image.
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

  if [[ "$path" == /work/* || "$path" == /usr/src/app/* || "$path" == /spliceai-cache/* ]]; then
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
    --precomputed-snv-vcf|--precomputed-indel-vcf|--precomputed-vcf|--genome)
      mapped_args+=("$1")
      if [[ $# -lt 2 ]]; then
        echo "Error: missing value for $1" >&2
        exit 1
      fi
      mapped_args+=("$(map_to_container_path "$2")")
      shift 2
      ;;
    --precomputed-snv-vcf=*|--precomputed-indel-vcf=*|--precomputed-vcf=*|--genome=*)
      key="${1%%=*}"
      val="${1#*=}"
      mapped_args+=("${key}=$(map_to_container_path "$val")")
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
  annotate-spliceai "$input_in_container" "$output_in_container" "${mapped_args[@]}"