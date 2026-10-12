#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 2 ]]; then
  cat <<'EOF'
Usage: src/scripts/run_liftover.sh <input.tsv> <output.tsv> [liftover options...]

Lift genomic positions between genome assemblies (e.g. hg19 -> hg38) using a
UCSC chain file. Two modes are available and may be combined:

  --input-col COL   HGVS mode: lift genomic HGVS g. strings in COL
  --chrom-col COL   Coordinate mode: lift chrom/start[/end] columns

Examples:
  src/scripts/run_liftover.sh variants.tsv lifted.tsv \
    --from-assembly hg19 --to-assembly hg38 --input-col hgvs_g

  src/scripts/run_liftover.sh regions.tsv lifted.tsv \
    --from-assembly hg19 --to-assembly hg38 \
    --chrom-col chrom --start-col start --end-col end

Notes:
  - Requires the pyliftover package: pip install variant-annotation[liftover]
  - Chain files are downloaded from UCSC on first use and cached in
    ~/.pyliftover/; pass --chain-file to supply a pre-downloaded file.
  - Runs locally using .venv/bin/python (not Docker).
  - See docs/liftover.md for full reference documentation.
EOF
  exit 1
fi

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/../.." && pwd)"

INPUT_TSV="$1"
OUTPUT_TSV="$2"
shift 2

cd "$REPO_DIR"

.venv/bin/python -m src.liftover "$INPUT_TSV" "$OUTPUT_TSV" "$@"
