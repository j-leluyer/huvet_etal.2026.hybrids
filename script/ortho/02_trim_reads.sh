#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 SAMPLE_BASE [PROJECT_DIR]" >&2
  exit 1
fi

SAMPLE_BASE="$1"
PROJECT_DIR="${2:-$(cd "$(dirname "$0")/../.." && pwd)}"
cd "$PROJECT_DIR"

command -v trimmomatic >/dev/null 2>&1 || { echo "trimmomatic is required in PATH" >&2; exit 1; }

RAW_DIR="${RAW_DIR:-raw_fastq}"
OUT_DIR="${OUT_DIR:-3-trimmed}"
ADAPTER_FILE="${ADAPTER_FILE:-adapters/univec.fasta}"
THREADS="${THREADS:-8}"

mkdir -p "$OUT_DIR"

R1="$RAW_DIR/${SAMPLE_BASE}_R1.fastq.gz"
R2="$RAW_DIR/${SAMPLE_BASE}_R2.fastq.gz"
[[ -f "$R1" ]] || { echo "Missing $R1" >&2; exit 1; }
[[ -f "$R2" ]] || { echo "Missing $R2" >&2; exit 1; }
[[ -f "$ADAPTER_FILE" ]] || { echo "Missing adapter file $ADAPTER_FILE" >&2; exit 1; }

trimmomatic PE -threads "$THREADS" -phred33 \
  "$R1" "$R2" \
  "$OUT_DIR/${SAMPLE_BASE}_R1.paired.fastq.gz" \
  "$OUT_DIR/${SAMPLE_BASE}_R1.single.fastq.gz" \
  "$OUT_DIR/${SAMPLE_BASE}_R2.paired.fastq.gz" \
  "$OUT_DIR/${SAMPLE_BASE}_R2.single.fastq.gz" \
  ILLUMINACLIP:"$ADAPTER_FILE":2:20:7 \
  LEADING:20 TRAILING:20 SLIDINGWINDOW:30:30 MINLEN:60
