#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  echo "Usage: $0 <input.vcf.gz|input.bcf> [out_dir]" >&2
  exit 1
fi

INPUT_VCF="$1"
OUT_DIR="${2:-output/vcf_ase_qc}"

command -v bcftools >/dev/null 2>&1 || { echo "bcftools is required in PATH" >&2; exit 1; }
command -v awk >/dev/null 2>&1 || { echo "awk is required in PATH" >&2; exit 1; }
command -v gzip >/dev/null 2>&1 || { echo "gzip is required in PATH" >&2; exit 1; }

mkdir -p "$OUT_DIR"

RAW_METRICS="$OUT_DIR/.vcf_raw_metrics.tsv"
SITE_METRICS="$OUT_DIR/vcf_site_metrics.tsv"
SUMMARY_OUT="$OUT_DIR/vcf_summary_stats.tsv"
GRID_OUT="$OUT_DIR/vcf_threshold_grid.tsv"
RECO_OUT="$OUT_DIR/vcf_recommended_thresholds.tsv"

# 1) Extract key per-site values; fill AN/AC/F_MISSING/DP tags when missing
bcftools +fill-tags "$INPUT_VCF" -Ou -- -t AN,AC,NS,F_MISSING,DP \
  | bcftools query -f '%CHROM\t%POS\t%QUAL\t%INFO/DP\t%AN\t%AC\t%INFO/F_MISSING\n' \
  > "$RAW_METRICS"

if [[ ! -s "$RAW_METRICS" ]]; then
  echo "No variants found in $INPUT_VCF" >&2
  exit 1
fi

# 2) Estimate sample count from max AN/2
N_SAMPLES_EST=$(awk 'BEGIN{m=0} $5!="." && $5!="" {v=$5+0; if(v>m) m=v} END{if(m>0) printf("%.0f", m/2); else print "NA"}' "$RAW_METRICS")

# 3) Build per-site metrics table
awk -v n_samples_est="$N_SAMPLES_EST" 'BEGIN{
  OFS="\t";
  print "CHROM","POS","QUAL","INFO_DP","AN","AC","AC_SUM","AF","F_MISSING","MISSING","EXP_HET";
}
{
  chrom=$1; pos=$2; qual=$3; dp=$4; an=$5; ac=$6; fmiss=$7;

  acsum="NA";
  if(ac != "." && ac != "" && ac != "NA") {
    n=split(ac, arr, ",");
    s=0; ok=0;
    for(i=1;i<=n;i++) {
      if(arr[i] != "." && arr[i] != "" && arr[i] != "NA") { s += arr[i]+0; ok=1; }
    }
    if(ok) acsum=s;
  }

  missing="NA";
  if(fmiss != "." && fmiss != "" && fmiss != "NA") {
    missing=fmiss+0;
  } else if(n_samples_est != "NA" && an != "." && an != "" && an != "NA") {
    missing = 1 - ((an+0) / (2*n_samples_est));
  }

  af="NA";
  if(acsum != "NA" && an != "." && an != "" && an != "NA" && (an+0)>0) {
    af = acsum/(an+0);
  }

  exphet="NA";
  if(af != "NA") exphet = 2*af*(1-af);

  print chrom,pos,qual,dp,an,ac,acsum,af,fmiss,missing,exphet;
}' "$RAW_METRICS" > "$SITE_METRICS"

gzip -f "$SITE_METRICS"
SITE_METRICS_GZ="$SITE_METRICS.gz"

# 4) Global summary (quick stats)
awk 'BEGIN{OFS="\t"}
NR==1{next}
{
  n++;

  if($4!="NA" && $4!="." && $4!="") {sDP+=$4; nDP++}
  if($10!="NA" && $10!="." && $10!="") {sMISS+=$10; nMISS++}
  if($5!="NA" && $5!="." && $5!="") {sAN+=$5; nAN++}
  if($7!="NA" && $7!="." && $7!="") {sAC+=$7; nAC++}
  if($8!="NA" && $8!="." && $8!="") {sAF+=$8; nAF++}
  if($11!="NA" && $11!="." && $11!="") {sHET+=$11; nHET++}
}
END{
  print "metric","value";
  print "n_variants",n;
  print "mean_DP", (nDP? sDP/nDP : "NA");
  print "mean_missing", (nMISS? sMISS/nMISS : "NA");
  print "mean_AN", (nAN? sAN/nAN : "NA");
  print "mean_AC_SUM", (nAC? sAC/nAC : "NA");
  print "mean_AF", (nAF? sAF/nAF : "NA");
  print "mean_expected_het", (nHET? sHET/nHET : "NA");
}' <(gzip -dc "$SITE_METRICS_GZ") > "$SUMMARY_OUT"

# 5) Threshold grid for ASE-oriented filtering
DP_LIST=(5 8 10 12 15)
MISS_LIST=(0.50 0.40 0.30 0.20 0.10)
AC_LIST=(1 2 3 4)

TOTAL=$(awk 'NR>1{n++} END{print n+0}' <(gzip -dc "$SITE_METRICS_GZ"))

{
  echo -e "dp_min\tmissing_max\tac_min\tn_retained\tprop_retained\tmean_expected_het\tinformative_score"
  for dp in "${DP_LIST[@]}"; do
    for miss in "${MISS_LIST[@]}"; do
      for acmin in "${AC_LIST[@]}"; do
        awk -v dp="$dp" -v miss="$miss" -v acmin="$acmin" -v total="$TOTAL" 'BEGIN{OFS="\t"}
        NR==1{next}
        {
          ok=1
          if($4=="NA" || $4=="." || $4=="" || ($4+0)<dp) ok=0
          if($10=="NA" || $10=="." || $10=="" || ($10+0)>miss) ok=0
          if($7=="NA" || $7=="." || $7=="" || ($7+0)<acmin) ok=0

          if(ok){
            n++
            if($11!="NA" && $11!="." && $11!="") shet += ($11+0)
          }
        }
        END{
          prop = (total>0 ? n/total : 0)
          mhet = (n>0 ? shet/n : "NA")
          score = (n>0 && mhet!="NA" ? n*mhet : "NA")
          print dp, miss, acmin, n+0, prop, mhet, score
        }' <(gzip -dc "$SITE_METRICS_GZ")
      done
    done
  done
} > "$GRID_OUT"

# sort by informative_score then n_retained (descending), keep header first
{
  head -n 1 "$GRID_OUT"
  tail -n +2 "$GRID_OUT" | sort -t $'\t' -k7,7gr -k4,4gr
} > "$GRID_OUT.tmp" && mv "$GRID_OUT.tmp" "$GRID_OUT"

# 6) Recommended profiles from grid
awk 'BEGIN{OFS="\t"; print "profile","dp_min","missing_max","ac_min","n_retained","prop_retained","mean_expected_het","informative_score"}
NR==1{next}
{
  dp=$1+0; miss=$2+0; ac=$3+0
  if(!p && miss<=0.50 && dp>=5  && ac>=1){ print "permissive",$0; p=1 }
  if(!b && miss<=0.30 && dp>=8  && ac>=2){ print "balanced",$0; b=1 }
  if(!s && miss<=0.20 && dp>=12 && ac>=3){ print "strict",$0; s=1 }
}
END{
  if(!p) print "permissive","NA","NA","NA","NA","NA","NA","NA"
  if(!b) print "balanced","NA","NA","NA","NA","NA","NA","NA"
  if(!s) print "strict","NA","NA","NA","NA","NA","NA","NA"
}' "$GRID_OUT" > "$RECO_OUT"

echo "Done. Files written to: $OUT_DIR"
echo " - $SITE_METRICS_GZ"
echo " - $SUMMARY_OUT"
echo " - $GRID_OUT"
echo " - $RECO_OUT"
