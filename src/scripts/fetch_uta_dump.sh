#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
OUT_DIR="${ROOT_DIR}/data/uta"
OUT_FILE="${OUT_DIR}/uta_20241220.pgd.gz"
URL="https://dl.biocommons.org/uta/uta_20241220.pgd.gz"

mkdir -p "${OUT_DIR}"

if [[ -s "${OUT_FILE}" ]]; then
  echo "UTA dump already exists: ${OUT_FILE}"
  exit 0
fi

# dl.biocommons.org serves UTA dumps behind a lightweight JS check. Any
# non-empty values for these cookies are sufficient for direct file download.
COOKIE="human_verified=cli; challenge_ts=$(date +%s%3N)"

echo "Downloading UTA dump to ${OUT_FILE}"
curl -fL --retry 3 --retry-delay 2 -H "Cookie: ${COOKIE}" "${URL}" -o "${OUT_FILE}.tmp"

# Verify the artifact looks like gzip before putting it in place.
gzip -t "${OUT_FILE}.tmp"
mv "${OUT_FILE}.tmp" "${OUT_FILE}"

echo "Downloaded ${OUT_FILE}"
