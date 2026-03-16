#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="${1:-$(cd "$(dirname "$0")/../.." && pwd)}"
shift || true
cd "$PROJECT_DIR"

MATRIX_SCRIPT=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --matrix-script)
      MATRIX_SCRIPT="$2"
      shift 2
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

for cmd in seqkit diamond TransDecoder.LongOrfs TransDecoder.Predict python3; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "$cmd is required in PATH" >&2; exit 1; }
done

TX2GENE="0_ref/combined_tx/tx2gene.with_mito.tsv"
CA_REF="0_ref/CA.ref.fa"
CG_REF="0_ref/CG.ref.fa"
ORTHO_IN="ortholog_1to1.tsv"
ORTHO_OUT="ortholog_1to1.with_mito.tsv"
CA_MITO_TX2GENE="0_ref/CA_mito_tx2gene.tsv"
CG_MITO_TX2GENE="0_ref/CG_mito_tx2gene.tsv"

[[ -f "$TX2GENE" ]] || { echo "Missing $TX2GENE" >&2; exit 1; }
[[ -f "$CA_REF" ]] || { echo "Missing $CA_REF" >&2; exit 1; }
[[ -f "$CG_REF" ]] || { echo "Missing $CG_REF" >&2; exit 1; }

awk '{print $1"\t"$2}' "$TX2GENE" | grep '^CA' > CA.tx2gene
seqkit fx2tab -n -l "$CA_REF" | awk -F'\t' 'BEGIN{OFS="\t"}{split($1,a," "); print a[1], $2}' > CA.tx.len
sort -k1,1 CA.tx2gene > CA.tx2gene.sorted
sort -k1,1 CA.tx.len > CA.tx.len.sorted
join -t $'\t' CA.tx2gene.sorted CA.tx.len.sorted | awk 'BEGIN{FS=OFS="\t"}{tx=$1; gene=$2; len=$3; if(!(gene in best) || len > bestlen[gene]){best[gene]=tx; bestlen[gene]=len}} END{for(g in best) print best[g], g}' > CA.rep_tx.tsv
cut -f1 CA.rep_tx.tsv > CA.rep_tx.ids
seqkit grep -f CA.rep_tx.ids "$CA_REF" > CA_rep_tx.fa

awk '{print $1"\t"$2}' "$TX2GENE" | grep '^CG' > CG.tx2gene
seqkit fx2tab -n -l "$CG_REF" | awk -F'\t' 'BEGIN{OFS="\t"}{split($1,a," "); print a[1], $2}' > CG.tx.len
sort -k1,1 CG.tx2gene > CG.tx2gene.sorted
sort -k1,1 CG.tx.len > CG.tx.len.sorted
join -t $'\t' CG.tx2gene.sorted CG.tx.len.sorted | awk 'BEGIN{FS=OFS="\t"}{tx=$1; gene=$2; len=$3; if(!(gene in best) || len > bestlen[gene]){best[gene]=tx; bestlen[gene]=len}} END{for(g in best) print best[g], g}' > CG.rep_tx.tsv
cut -f1 CG.rep_tx.tsv > CG.rep_tx.ids
seqkit grep -f CG.rep_tx.ids "$CG_REF" > CG_rep_tx.fa

TransDecoder.LongOrfs -t CG_rep_tx.fa
TransDecoder.Predict -t CG_rep_tx.fa
cp CG_rep_tx.fa.transdecoder.pep CG_bestORF_prot.fa

TransDecoder.LongOrfs -t CA_rep_tx.fa
TransDecoder.Predict -t CA_rep_tx.fa
cp CA_rep_tx.fa.transdecoder.pep CA_bestORF_prot.fa

diamond makedb --in CG_bestORF_prot.fa -d CG.db
diamond makedb --in CA_bestORF_prot.fa -d CA.db

diamond blastp -q CA_bestORF_prot.fa -d CG.db -o CA_vs_CG.tsv --max-target-seqs 25 --evalue 1e-10 -p 12 --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore
diamond blastp -q CG_bestORF_prot.fa -d CA.db -o CG_vs_CA.tsv --max-target-seqs 25 --evalue 1e-10 -p 12 --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore

LC_ALL=C sort -t $'\t' -k1,1 -k8,8nr -k7,7g CA_vs_CG.tsv | awk -F'\t' 'BEGIN{OFS="\t"} ($4/$5)>=0.5 && ($4/$6)>=0.5 && !seen[$1]++ {print}' > CA_vs_CG.best.tsv
LC_ALL=C sort -t $'\t' -k1,1 -k8,8nr -k7,7g CG_vs_CA.tsv | awk -F'\t' 'BEGIN{OFS="\t"} ($4/$5)>=0.5 && ($4/$6)>=0.5 && !seen[$1]++ {print}' > CG_vs_CA.best.tsv

cut -f1-2 CA_vs_CG.best.tsv | LC_ALL=C sort > A.best.pairs
awk -F'\t' 'BEGIN{OFS="\t"}{print $2,$1}' CG_vs_CA.best.tsv | LC_ALL=C sort > B.best.pairs_swapped
comm -12 A.best.pairs B.best.pairs_swapped > RBH_1to1.tsv

if [[ -n "$MATRIX_SCRIPT" ]]; then
  [[ -f "$MATRIX_SCRIPT" ]] || { echo "Missing matrix script $MATRIX_SCRIPT" >&2; exit 1; }
  python3 "$MATRIX_SCRIPT"
else
  echo "Skipping ortholog matrix helper: pass --matrix-script /path/to/make_ortholog_matrix.py to build $ORTHO_IN" >&2
fi

[[ -f "$ORTHO_IN" ]] || { echo "Missing $ORTHO_IN after ortholog build" >&2; exit 1; }
[[ -f "$CA_MITO_TX2GENE" ]] || { echo "Missing $CA_MITO_TX2GENE" >&2; exit 1; }
[[ -f "$CG_MITO_TX2GENE" ]] || { echo "Missing $CG_MITO_TX2GENE" >&2; exit 1; }

awk -F'\t' 'NR==1{next} $2 ~ /^CAmt\|/ {split($2,a,"|"); print a[2]} $3 ~ /^CGmt\|/ {split($3,a,"|"); print a[2]}' "$ORTHO_IN" | LC_ALL=C sort -u > existing_mito.keys
awk -F'\t' '$2 ~ /^CAmt\|/ {split($2,a,"|"); print a[2]}' "$CA_MITO_TX2GENE" | LC_ALL=C sort -u > CAmt.keys
awk -F'\t' '$2 ~ /^CGmt\|/ {split($2,a,"|"); print a[2]}' "$CG_MITO_TX2GENE" | LC_ALL=C sort -u > CGmt.keys
LC_ALL=C comm -12 CAmt.keys CGmt.keys > mito.keys.common
LC_ALL=C comm -23 mito.keys.common existing_mito.keys > mito.keys.to_add
awk 'BEGIN{OFS="\t"} {print "CAmt|"$1, "CGmt|"$1}' mito.keys.to_add > mito_gene_pairs.to_add.tsv
max_id_num=$(awk -F'\t' 'NR>1{gsub(/^ORTHO_/,"",$1); if($1+0>m)m=$1+0} END{print m+0}' "$ORTHO_IN")
start=$((max_id_num + 1))
awk -F'\t' -v OFS='\t' -v start="$start" '{ printf "ORTHO_%06d\t%s\t%s\n", start+NR-1, $1, $2 }' mito_gene_pairs.to_add.tsv > mito_ortholog_append.tsv
cat "$ORTHO_IN" mito_ortholog_append.tsv > "$ORTHO_OUT"

echo "Created $ORTHO_OUT"
