#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
cd "$PROJECT_DIR"

BASE="${BASE:-4-mapped-salmon}"
ORTHO="${ORTHO:-ortholog_1to1.with_mito.tsv}"
TX2GENE="${TX2GENE:-info/tx2gene.with_mito.tsv}"
OUTDIR="${OUTDIR:-ortho_quants}"
OUTMATRIX="${OUTMATRIX:-ortho_counts_matrix.tsv}"

[[ -f "$ORTHO" ]] || { echo "Missing $ORTHO" >&2; exit 1; }
[[ -f "$TX2GENE" ]] || { echo "Missing $TX2GENE" >&2; exit 1; }
[[ -d "$BASE" ]] || { echo "Missing $BASE" >&2; exit 1; }

mkdir -p "$OUTDIR"

awk 'BEGIN{FS=OFS="\t"} NR>1 {print $2,$1; print $3,$1}' "$ORTHO" | tr -d '\r' | awk 'NF==2' | LC_ALL=C sort -u > gene2ortho.tsv
LC_ALL=C sort -t $'\t' -k1,1 gene2ortho.tsv > gene2ortho.sorted.tsv

awk 'BEGIN{OFS="\t"} {gsub(/\r/,""); if(NF>=2) print $1,$2}' "$TX2GENE" | LC_ALL=C sort -t $'\t' -k2,2 > tx2gene.by_gene.tsv
join -t $'\t' -1 2 -2 1 tx2gene.by_gene.tsv gene2ortho.sorted.tsv | awk -F'\t' 'BEGIN{OFS="\t"} {print $2,$3}' | LC_ALL=C sort -u > tx2ortho.tsv

for d in "$BASE"/*; do
  [[ -d "$d" ]] || continue
  s=$(basename "$d")
  q="$d/quant.sf"
  [[ -f "$q" ]] || continue

  gawk -v FS='\t' -v OFS='\t' '
    FNR==NR { tx2o[$1]=$2; next }
    FNR==1 { next }
    {
      o = tx2o[$1]
      if (o != "") {
        tpm[o] += $4
        reads[o] += $5
      }
    }
    END { for (o in reads) print o, tpm[o], reads[o] }
  ' tx2ortho.tsv "$q" | LC_ALL=C sort -t $'\t' -k1,1 > "$OUTDIR/${s}.ortho.tsv"
done

awk 'BEGIN{FS="\t"} NR>1{print $1}' "$ORTHO" | tr -d '\r' | LC_ALL=C sort -u > ORTHO_master.list

{
  printf "ORTHO_ID"
  for f in "$OUTDIR"/*.ortho.tsv; do
    printf "\t%s" "$(basename "$f" .ortho.tsv)"
  done
  printf "\n"
} > "$OUTMATRIX"

gawk 'BEGIN{FS=OFS="\t"}
  FNR==NR { master[++n]=$1; next }
  FNR==1 { ns++ }
  { data[$1,ns]=$3 }
  END{
    for(i=1;i<=n;i++){
      o=master[i]
      printf "%s", o
      for(j=1;j<=ns;j++){
        v=data[o,j]
        if(v=="") v=0
        printf OFS v
      }
      printf "\n"
    }
  }' ORTHO_master.list "$OUTDIR"/*.ortho.tsv >> "$OUTMATRIX"

echo "Created $OUTMATRIX"
