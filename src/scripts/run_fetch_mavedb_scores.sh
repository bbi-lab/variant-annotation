#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  cat <<'EOF'
Usage: src/scripts/run_fetch_mavedb_scores.sh <output-file> [fetch_mavedb_scores options...]

Examples:
  src/scripts/run_fetch_mavedb_scores.sh data/excalibr_domainome/excalibr_domainome.0.tsv \
    --score-sets-file data/excalibr_domainome/score_sets.tsv
  src/scripts/run_fetch_mavedb_scores.sh /work/variants.0.tsv \
    --score-set-urn urn:mavedb:00000895-a-1 \
    --score-set-urn urn:mavedb:00000736-a-1

Notes:
  - A relative output path is written under the repository (/usr/src/app in the
    container) if its parent directory exists on the host, otherwise under /work.
  - By default /work maps to ./data on the host.
  - Override mount directory with VARIANT_DATA_DIR=/absolute/path.
  - Option paths such as --score-sets-file are resolved relative to the
    repository root.
  - The MaveDB API URL defaults to https://api.mavedb.org.  Set MAVEDB_API_URL
    in settings/.env to override (e.g. to point at a local instance).
  - Images are reused by default for fast runs.
  - Add --rebuild-image to force rebuilding the image.
  - Add --no-build-cache with --rebuild-image for a clean rebuild.

Output columns:
  gene_symbol, score_set_title, variant_urn, raw_hgvs_nt, raw_hgvs_pro,
  score, rna_score, rna_score_d6, rna_score_d20, score_data, reference_based,
  target_sequence_type, target_sequence, preferred_transcript
  With --include-mavedb-mappings, also mavedb_mapped_hgvs_g,
  mavedb_mapped_hgvs_c, and mavedb_mapped_hgvs_p (after reference_based).
EOF
  exit 1
fi

output_path="$1"
shift 1

compose_build_flag=""
compose_no_cache_flag=""

# Wrapper-only options (must appear before fetch_mavedb_scores options).
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

if [[ "$output_path" == /work/* || "$output_path" == /usr/src/app/* ]]; then
  output_in_container="$output_path"
elif [[ -d "$(dirname "$output_path")" ]]; then
  output_in_container="/usr/src/app/$output_path"
else
  output_in_container="/work/$output_path"
fi

exec docker compose --profile tools run \
  ${compose_build_flag:+$compose_build_flag} ${compose_no_cache_flag:+$compose_no_cache_flag} --rm \
  fetch-mavedb-scores "$output_in_container" "$@"
