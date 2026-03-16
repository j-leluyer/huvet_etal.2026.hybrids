#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 SAMPLE_BASE [PROJECT_DIR]" >&2
  exit 1
fi

SAMPLE_BASE="$1"
PROJECT_DIR="${2:-$(cd "$(dirname "$0")/../.." && pwd)}"
cd "$PROJECT_DIR"

command -v salmon >/dev/null 2>&1 || { echo "salmon is required in PATH" >&2; exit 1; }

INDEX_DIR="${INDEX_DIR:-0_ref/combined_tx/salmon_index}"
TRIM_DIR="${TRIM_DIR:-3-trimmed}"
OUT_DIR="${OUT_DIR:-4-mapped-salmon}"
THREADS="${THREADS:-12}"

R1="$TRIM_DIR/${SAMPLE_BASE}_R1.paired.fastq.gz"
R2="$TRIM_DIR/${SAMPLE_BASE}_R2.paired.fastq.gz"
[[ -d "$INDEX_DIR" ]] || { echo "Missing salmon index $INDEX_DIR" >&2; exit 1; }
[[ -f "$R1" ]] || { echo "Missing $R1" >&2; exit 1; }
[[ -f "$R2" ]] || { echo "Missing $R2" >&2; exit 1; }
mkdir -p "$OUT_DIR"

salmon quant \
  -i "$INDEX_DIR" \
  -l A \
  -1 "$R1" \
  -2 "$R2" \
  --validateMappings \
  --seqBias \
  --gcBias \
  -p "$THREADS" \
  -o "$OUT_DIR/$SAMPLE_BASE"
