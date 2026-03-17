#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  cat(
    "Usage:\n",
    "  Rscript script/summarize_vcf_ase_thresholds.R <input.vcf.gz|input.bcf> [out_dir]\n\n",
    "Outputs:\n",
    "  - vcf_site_metrics.tsv.gz\n",
    "  - vcf_summary_stats.tsv\n",
    "  - vcf_threshold_grid.tsv\n",
    "  - vcf_recommended_thresholds.tsv\n",
    sep = ""
  )
  quit(status = 1)
}

input_vcf <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- if (length(args) >= 2) args[[2]] else file.path("output", "vcf_ase_qc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

if (Sys.which("bcftools") == "") stop("bcftools not found in PATH")

tmp_tsv <- tempfile(pattern = "vcf_metrics_", fileext = ".tsv")
on.exit(unlink(tmp_tsv), add = TRUE)

# Fill key INFO tags to avoid missing AC/AN/F_MISSING when absent in input VCF
query_fmt <- "%CHROM\\t%POS\\t%QUAL\\t%INFO/DP\\t%AN\\t%AC\\t%INFO/F_MISSING\\n"
cmd <- paste(
  "bcftools +fill-tags", shQuote(input_vcf), "-Ou -- -t AN,AC,NS,F_MISSING,DP",
  "| bcftools query -f", shQuote(query_fmt), ">", shQuote(tmp_tsv)
)
status <- system(cmd)
if (status != 0) stop("Failed to extract metrics from VCF with bcftools")

x <- fread(
  tmp_tsv,
  header = FALSE,
  sep = "\t",
  na.strings = c(".", "NA", ""),
  col.names = c("CHROM", "POS", "QUAL", "INFO_DP", "AN", "AC", "F_MISSING")
)

if (nrow(x) == 0) stop("No variants found in input VCF")

sum_ac <- function(v) {
  if (is.na(v) || v == "") return(NA_real_)
  vals <- suppressWarnings(as.numeric(strsplit(v, ",", fixed = TRUE)[[1]]))
  if (all(is.na(vals))) return(NA_real_)
  sum(vals, na.rm = TRUE)
}

x[, AC_SUM := vapply(AC, sum_ac, numeric(1))]
x[, QUAL := as.numeric(QUAL)]
x[, INFO_DP := as.numeric(INFO_DP)]
x[, AN := as.numeric(AN)]
x[, F_MISSING := as.numeric(F_MISSING)]

n_samples_est <- suppressWarnings(max(x$AN, na.rm = TRUE) / 2)
if (!is.finite(n_samples_est) || is.na(n_samples_est) || n_samples_est <= 0) {
  n_samples_est <- NA_real_
}

x[, MISSING_FROM_AN := if (is.finite(n_samples_est)) 1 - (AN / (2 * n_samples_est)) else NA_real_]
x[, MISSING := fifelse(!is.na(F_MISSING), F_MISSING, MISSING_FROM_AN)]
x[, AF := fifelse(!is.na(AN) & AN > 0 & !is.na(AC_SUM), AC_SUM / AN, NA_real_)]
x[, EXP_HET := fifelse(!is.na(AF), 2 * AF * (1 - AF), NA_real_)]

# Save per-site metrics (compressed)
site_out <- file.path(out_dir, "vcf_site_metrics.tsv.gz")
fwrite(
  x[, .(CHROM, POS, QUAL, INFO_DP, AN, AC, AC_SUM, AF, F_MISSING, MISSING, EXP_HET)],
  file = site_out,
  sep = "\t"
)

qf <- function(v, p) as.numeric(quantile(v, probs = p, na.rm = TRUE, names = FALSE, type = 7))

summary_tbl <- data.table(
  n_variants = nrow(x),
  n_samples_est = n_samples_est,
  mean_DP = mean(x$INFO_DP, na.rm = TRUE),
  median_DP = median(x$INFO_DP, na.rm = TRUE),
  dp_q10 = qf(x$INFO_DP, 0.10),
  dp_q25 = qf(x$INFO_DP, 0.25),
  dp_q50 = qf(x$INFO_DP, 0.50),
  dp_q75 = qf(x$INFO_DP, 0.75),
  dp_q90 = qf(x$INFO_DP, 0.90),
  mean_missing = mean(x$MISSING, na.rm = TRUE),
  missing_q50 = qf(x$MISSING, 0.50),
  missing_q75 = qf(x$MISSING, 0.75),
  missing_q90 = qf(x$MISSING, 0.90),
  mean_AN = mean(x$AN, na.rm = TRUE),
  median_AN = median(x$AN, na.rm = TRUE),
  mean_AC_SUM = mean(x$AC_SUM, na.rm = TRUE),
  median_AC_SUM = median(x$AC_SUM, na.rm = TRUE),
  mean_AF = mean(x$AF, na.rm = TRUE),
  median_AF = median(x$AF, na.rm = TRUE),
  mean_expected_het = mean(x$EXP_HET, na.rm = TRUE)
)

fwrite(summary_tbl, file.path(out_dir, "vcf_summary_stats.tsv"), sep = "\t")

# Threshold grid for ASE-oriented retention/informativeness trade-off
DP_MIN <- c(5, 8, 10, 12, 15)
MISS_MAX <- c(0.50, 0.40, 0.30, 0.20, 0.10)
AC_MIN <- c(1, 2, 3, 4)

grid <- CJ(dp_min = DP_MIN, missing_max = MISS_MAX, ac_min = AC_MIN)

n_total <- nrow(x)
grid[, c("n_retained", "prop_retained", "mean_expected_het", "informative_score") := {
  keep <- !is.na(INFO_DP) & INFO_DP >= dp_min &
    !is.na(MISSING) & MISSING <= missing_max &
    !is.na(AC_SUM) & AC_SUM >= ac_min

  n_keep <- sum(keep)
  prop <- n_keep / n_total
  mhet <- if (n_keep > 0) mean(x$EXP_HET[keep], na.rm = TRUE) else NA_real_
  score <- if (n_keep > 0 && is.finite(mhet)) n_keep * mhet else NA_real_

  list(n_keep, prop, mhet, score)
}, by = .(dp_min, missing_max, ac_min)]

setorder(grid, -informative_score, -n_retained)
fwrite(grid, file.path(out_dir, "vcf_threshold_grid.tsv"), sep = "\t")

pick_best <- function(dt, cond) {
  sub <- dt[eval(cond)]
  if (nrow(sub) == 0) return(data.table(profile = NA_character_, dp_min = NA_real_, missing_max = NA_real_, ac_min = NA_real_, n_retained = NA_real_, prop_retained = NA_real_, mean_expected_het = NA_real_, informative_score = NA_real_))
  sub[order(-sub$informative_score, -sub$n_retained)][1]
}

best_permissive <- pick_best(grid, quote(missing_max <= 0.50 & dp_min >= 5 & ac_min >= 1))
best_balanced <- pick_best(grid, quote(missing_max <= 0.30 & dp_min >= 8 & ac_min >= 2))
best_strict <- pick_best(grid, quote(missing_max <= 0.20 & dp_min >= 12 & ac_min >= 3))

best_permissive[, profile := "permissive"]
best_balanced[, profile := "balanced"]
best_strict[, profile := "strict"]

reco <- rbindlist(list(best_permissive, best_balanced, best_strict), fill = TRUE)
setcolorder(reco, c("profile", "dp_min", "missing_max", "ac_min", "n_retained", "prop_retained", "mean_expected_het", "informative_score"))
fwrite(reco, file.path(out_dir, "vcf_recommended_thresholds.tsv"), sep = "\t")

cat("Done. Files written to:", out_dir, "\n")
cat(" -", file.path(out_dir, "vcf_site_metrics.tsv.gz"), "\n")
cat(" -", file.path(out_dir, "vcf_summary_stats.tsv"), "\n")
cat(" -", file.path(out_dir, "vcf_threshold_grid.tsv"), "\n")
cat(" -", file.path(out_dir, "vcf_recommended_thresholds.tsv"), "\n")
