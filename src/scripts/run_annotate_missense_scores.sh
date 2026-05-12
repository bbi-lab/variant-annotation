#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  cat <<'EOF'
Usage: src/scripts/run_annotate_missense_scores.sh <input-file> <output-file> [annotate_missense_scores options...]

Examples:
  src/scripts/run_annotate_missense_scores.sh variants.tsv annotated.tsv \
    --revel-file /path/to/revel_hg38.tsv.gz \
    --alphamissense-file /path/to/AlphaMissense_hg38.tsv.gz \
    --dbnsfp-file /path/to/dbNSFP5.3.1a_grch38.gz

  # Using environment variables for file paths:
  export REVEL_FILE=/path/to/revel_hg38.tsv.gz
  export ALPHAMISSENSE_FILE=/path/to/AlphaMissense_hg38.tsv.gz
  export DBNSFP_FILE=/path/to/dbNSFP5.3.1a_grch38.gz
  src/scripts/run_annotate_missense_scores.sh variants.tsv annotated.tsv

Notes:
  - Input/output paths are interpreted relative to /work in the container.
  - By default /work maps to ./data on the host.
  - Override mount directory with VARIANT_DATA_DIR=/absolute/path.
  - Score files (--revel-file, --alphamissense-file, --dbnsfp-file) are passed
    directly; if they are inside ./data they will be accessible under /work.
  - Requires tabix (htslib) in the container image.
  - At least one of --revel-file, --alphamissense-file, or --dbnsfp-file must
    be provided.
  - Images are reused by default for fast runs.
  - Add --rebuild-image to force rebuilding the image.
  - Add --no-build-cache with --rebuild-image for a clean rebuild.

Data file preparation:
  AlphaMissense (index must be generated locally):
    wget https://storage.googleapis.com/dm_alphamissense/AlphaMissense_hg38.tsv.gz
    tabix -s 1 -b 2 -e 2 AlphaMissense_hg38.tsv.gz

  dbNSFP (GRCh38 variant file; .tbi is available pre-built):
    wget https://dist.genos.us/academic/e55b09/dbNSFP5.3.1a_grch38.gz
    wget https://dist.genos.us/academic/e55b09/dbNSFP5.3.1a_grch38.gz.tbi

  REVEL (requires one-time preparation):
    wget https://rothsj06.dmz.hpc.mssm.edu/revel-v1.3_all_chromosomes.zip
    unzip revel-v1.3_all_chromosomes.zip
    mv revel_with_transcript_ids revel_with_transcript_ids.csv
    tail -n +2 revel_with_transcript_ids.csv \
      | awk -F',' 'NF>=9 && $3!="" && $3!="." {print $1"\t"$3"\t"$4"\t"$5"\t"$8}' \
      | (printf '#chr\tpos\tref\talt\trevel_score\n'; sort -k1,1V -k2,2n) \
      | bgzip > revel_hg38.tsv.gz
    tabix -s 1 -b 2 -e 2 -S 1 revel_hg38.tsv.gz
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
  annotate-missense-scores "$input_in_container" "$output_in_container" "$@"
