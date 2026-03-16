#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
cd "$PROJECT_DIR"

command -v salmon >/dev/null 2>&1 || { echo "salmon is required in PATH" >&2; exit 1; }

REFERENCE="${REFERENCE:-0_ref/combined_tx/CA_CG_all.tag.fa}"
INDEX_DIR="${INDEX_DIR:-0_ref/combined_tx/salmon_index}"
THREADS="${THREADS:-4}"
KMER="${KMER:-31}"

[[ -f "$REFERENCE" ]] || { echo "Missing reference FASTA $REFERENCE" >&2; exit 1; }

salmon index -t "$REFERENCE" -i "$INDEX_DIR" -k "$KMER" -p "$THREADS"
