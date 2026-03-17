#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 OUTPUT_FILE [PROJECT_DIR]" >&2
  exit 1
fi

OUT_FILE="$1"
PROJECT_DIR="${2:-$(cd "$(dirname "$0")/../.." && pwd)}"
cd "$PROJECT_DIR"

mkdir -p "$(dirname "$OUT_FILE")"

{
  echo "# Tool version snapshot"
  echo "date: $(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  echo

  log_cmd_version() {
    local label="$1"
    shift
    if command -v "$1" >/dev/null 2>&1; then
      local version
      version="$($@ 2>&1 | head -n 1 || true)"
      echo "$label: ${version:-available (version string not detected)}"
    else
      echo "$label: NOT_FOUND"
    fi
  }

  log_cmd_version "bedtools" bedtools --version
  log_cmd_version "trimmomatic" trimmomatic -version
  log_cmd_version "salmon" salmon --version
  log_cmd_version "seqkit" seqkit version
  log_cmd_version "diamond" diamond version
  log_cmd_version "python3" python3 --version

  if command -v TransDecoder.LongOrfs >/dev/null 2>&1; then
    td_ver="$(TransDecoder.LongOrfs --version 2>&1 | head -n 1 || true)"
    echo "TransDecoder.LongOrfs: ${td_ver:-available}"
  else
    echo "TransDecoder.LongOrfs: NOT_FOUND"
  fi

  if command -v TransDecoder.Predict >/dev/null 2>&1; then
    tdp_ver="$(TransDecoder.Predict --version 2>&1 | head -n 1 || true)"
    echo "TransDecoder.Predict: ${tdp_ver:-available}"
  else
    echo "TransDecoder.Predict: NOT_FOUND"
  fi
} > "$OUT_FILE"

echo "Wrote tool versions to $OUT_FILE"
