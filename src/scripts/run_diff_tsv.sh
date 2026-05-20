#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  cat <<'EOF'
Usage: src/scripts/run_diff_tsv.sh <old.tsv> <new.tsv> [output_prefix]

Row-level diff between two TSV files joined on variant_urn. Produces two
output files:
  <prefix>.changed_rows.tsv  — rows from the new file that are new or changed
  <prefix>.hgvs_diffs.tsv    — per-component diffs of mapped_hgvs_g/c/p

Arguments:
  old.tsv        Baseline (old) TSV file
  new.tsv        Updated (new) TSV file
  output_prefix  Optional prefix for output file names
                 Default: stem of new.tsv in the current directory

Notes:
  - Both files must contain a variant_urn column.
  - Only rows present in the new file are reported; rows only in the old
    file are not included in the output.
  - Runs locally using .venv/bin/python (not Docker).
EOF
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

OLD_TSV="$1"
NEW_TSV="$2"
OUTPUT_PREFIX="${3:-}"

cd "$REPO_DIR"

ARGS=("$OLD_TSV" "$NEW_TSV")
if [[ -n "$OUTPUT_PREFIX" ]]; then
    ARGS+=(-o "$OUTPUT_PREFIX")
fi

.venv/bin/python src/diff_tsv.py "${ARGS[@]}"
