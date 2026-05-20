#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  cat <<'EOF'
Usage: src/scripts/run_diff_hgvs.sh <old.tsv> <new.tsv> [output.tsv]

Compare the mapped_hgvs_g, mapped_hgvs_c, and mapped_hgvs_p columns between
two TSV files joined on variant_urn. Reports one output row per changed
pipe-separated component.

Arguments:
  old.tsv     Baseline (old) TSV file
  new.tsv     Updated (new) TSV file
  output.tsv  Optional output path
              Default: <new_tsv_stem>.hgvs_diffs.tsv in the current directory

Notes:
  - Both files must contain a variant_urn column.
  - Runs locally using .venv/bin/python (not Docker).
EOF
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

OLD_TSV="$1"
NEW_TSV="$2"
OUTPUT="${3:-}"

cd "$REPO_DIR"

ARGS=("$OLD_TSV" "$NEW_TSV")
if [[ -n "$OUTPUT" ]]; then
    ARGS+=(-o "$OUTPUT")
fi

.venv/bin/python src/diff_hgvs.py "${ARGS[@]}"
