#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
cd "$PROJECT_DIR"

command -v bedtools >/dev/null 2>&1 || { echo "bedtools is required in PATH" >&2; exit 1; }

CA_GFF="0_ref/genomic.CA.fixed.gff"
CG_GFF="0_ref/genomic.CG.fixed.gff"
COMBINED_DIR="0_ref/combined_tx"
COMBINED_FASTA="$COMBINED_DIR/CA_CG_all.tag.fa"
TX2GENE_OUT="$COMBINED_DIR/tx2gene.with_mito.tsv"

mkdir -p "$COMBINED_DIR"
[[ -f "$COMBINED_FASTA" ]] || echo "Warning: $COMBINED_FASTA not found; downstream sanity checks may fail." >&2

awk -F'\t' '
BEGIN{OFS="\t"}
$0 ~ /^#/ {next}
($3=="exon" || $3=="CDS" || $3=="gene" || $3=="region") {next}
{
  n=split($9, kv, ";")
  delete a
  for (i=1; i<=n; i++) {
    split(kv[i], p, "=")
    a[p[1]] = p[2]
  }
  tx=""
  if (("transcript_id" in a) && a["transcript_id"] ~ /^(X[MR]_|N[MR]_)/) tx=a["transcript_id"]
  else if (("Name" in a) && a["Name"] ~ /^(X[MR]_|N[MR]_)/) tx=a["Name"]
  else if ("Dbxref" in a) {
    if (match(a["Dbxref"], /GenBank:(X[MR]_[^,;]+)/, m)) tx=m[1]
    else if (match(a["Dbxref"], /GenBank:(N[MR]_[^,;]+)/, m)) tx=m[1]
  }
  g=""
  if ("gene" in a) g=a["gene"]
  else if ("Parent" in a) g=a["Parent"]
  sub(/^gene-/, "", g)
  if (tx != "" && g != "") print "CA|"tx, g
}
' "$CA_GFF" | sort -u > "$COMBINED_DIR/CA_tx2gene.tsv"

awk -F'\t' '
BEGIN{OFS="\t"}
$0 ~ /^#/ {next}
($3=="exon" || $3=="CDS" || $3=="gene" || $3=="region") {next}
{
  n=split($9, kv, ";")
  delete a
  for (i=1; i<=n; i++) {
    split(kv[i], p, "=")
    a[p[1]] = p[2]
  }
  tx=""
  if (("transcript_id" in a) && a["transcript_id"] ~ /^(X[MR]_|N[MR]_)/) tx=a["transcript_id"]
  else if (("Name" in a) && a["Name"] ~ /^(X[MR]_|N[MR]_)/) tx=a["Name"]
  else if ("Dbxref" in a) {
    if (match(a["Dbxref"], /GenBank:(X[MR]_[^,;]+)/, m)) tx=m[1]
    else if (match(a["Dbxref"], /GenBank:(N[MR]_[^,;]+)/, m)) tx=m[1]
  }
  g=""
  if ("gene" in a) g=a["gene"]
  else if ("Parent" in a) g=a["Parent"]
  sub(/^gene-/, "", g)
  if (tx != "" && g != "") print "CG|"tx, g
}
' "$CG_GFF" | sort -u > "$COMBINED_DIR/CG_tx2gene.tsv"

cat "$COMBINED_DIR/CA_tx2gene.tsv" "$COMBINED_DIR/CG_tx2gene.tsv" > "$COMBINED_DIR/tx2gene.tsv"

if [[ -f "$COMBINED_FASTA" ]]; then
  grep '^>' "$COMBINED_FASTA" | sed 's/^>//' | awk '{print $1}' | tr -d '\r' | LC_ALL=C sort -u > "$COMBINED_DIR/all_tx.ids"
  cut -f1 "$COMBINED_DIR/tx2gene.tsv" | tr -d '\r' | LC_ALL=C sort -u > "$COMBINED_DIR/tx2gene.ids"
  LC_ALL=C comm -23 "$COMBINED_DIR/all_tx.ids" "$COMBINED_DIR/tx2gene.ids" > "$COMBINED_DIR/missing_tx.ids"
  grep -E '^(CAmt\||CGmt\|)' "$COMBINED_DIR/missing_tx.ids" > "$COMBINED_DIR/missing_mito.ids" || true
  awk 'BEGIN{OFS="\t"} {print $0, $0}' "$COMBINED_DIR/missing_mito.ids" > "$COMBINED_DIR/mito.tx2gene.add.tsv"
  cat "$COMBINED_DIR/tx2gene.tsv" "$COMBINED_DIR/mito.tx2gene.add.tsv" | LC_ALL=C sort -u > "$TX2GENE_OUT"
  echo "Created $TX2GENE_OUT"
else
  cp "$COMBINED_DIR/tx2gene.tsv" "$TX2GENE_OUT"
  echo "Created $TX2GENE_OUT without FASTA/mito sanity augmentation"
fi
