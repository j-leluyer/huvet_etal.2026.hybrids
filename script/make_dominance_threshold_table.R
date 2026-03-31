#!/usr/bin/env Rscript
# Generates Table_S2.dominance_thresholds.xlsx
# Describes the parameters and classification rules used by compute_dominance_effects.R

suppressPackageStartupMessages({
  if (!requireNamespace("writexl", quietly = TRUE)) install.packages("writexl", repos = "https://cloud.r-project.org")
  library(writexl)
})

args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", args[grep("--file=", args)])
script_dir <- if (length(script_path) > 0) {
  dirname(normalizePath(script_path[1], mustWork = TRUE))
} else {
  getwd()
}

project_dir <- if (length(commandArgs(trailingOnly = TRUE)) >= 1) {
  normalizePath(commandArgs(trailingOnly = TRUE)[1], mustWork = TRUE)
} else {
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
}

setwd(project_dir)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------
# Sheet 1: Classification parameters
# ------------------------------------------------------------------
parameters <- data.frame(
  parameter = c(
    "alpha",
    "delta"
  ),
  value = c(
    "0.05",
    "0.25"
  ),
  description = c(
    "Significance threshold applied to DESeq2 Benjamini-Hochberg adjusted p-value (padj). Both the hybrid-vs-AA and parent GG-vs-AA contrasts must satisfy padj < alpha to be assigned a non-Ambiguous dominance class.",
    "Maximum allowed deviation in log2 fold-change units used to define class boundaries. A gene is assigned to a class only if its log2FC falls within delta of the expected expression level for that class."
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------
# Sheet 2: Classification rules
# All log2FC values are relative to AA as reference (DESeq2 contrast).
# Hybrid contrast: AG vs AA, or GA vs AA.
# Parent contrast: GG vs AA.
# ------------------------------------------------------------------
rules <- data.frame(
  dominance_class = c(
    "Additive",
    "GG_dominant",
    "AA_dominant",
    "Overdominant",
    "Underdominant",
    "Ambiguous"
  ),
  condition = c(
    "hybrid padj < alpha  AND  parent (GG vs AA) padj < alpha  AND  |log2FC_hybrid - 0.5 * log2FC_GG_vs_AA| < delta",
    "hybrid padj < alpha  AND  |log2FC_hybrid - log2FC_GG_vs_AA| < delta",
    "hybrid padj < alpha  AND  |log2FC_hybrid| < delta",
    "hybrid padj < alpha  AND  parent (GG vs AA) padj < alpha  AND  log2FC_hybrid > log2FC_GG_vs_AA + delta",
    "hybrid padj < alpha  AND  parent (GG vs AA) padj < alpha  AND  log2FC_hybrid < -delta",
    "None of the above conditions met"
  ),
  biological_interpretation = c(
    "Hybrid expression is intermediate between the two parents (mid-parental value). The hybrid log2FC is within delta of the expected mid-parent value (0.5 x log2FC_GG_vs_AA).",
    "Hybrid expression resembles the GG parent. The hybrid log2FC is within delta of the GG-vs-AA log2FC.",
    "Hybrid expression resembles the AA parent (reference level). The hybrid log2FC is within delta of zero.",
    "Hybrid expression exceeds the high parent (GG). Hybrid log2FC is more than delta above the GG-vs-AA log2FC (transgressive up-regulation).",
    "Hybrid expression falls below the low parent (AA). Hybrid log2FC is below -delta (transgressive down-regulation).",
    "Gene does not meet significance or log2FC proximity criteria for any class."
  ),
  stringsAsFactors = FALSE
)

# ------------------------------------------------------------------
# Sheet 3: Parameter values used (explicit numeric summary)
# ------------------------------------------------------------------
summary_values <- data.frame(
  item = c(
    "alpha (padj threshold)",
    "delta (log2FC tolerance)",
    "Mid-parent expected log2FC",
    "Additive window",
    "GG-dominant window",
    "AA-dominant window",
    "Overdominant threshold",
    "Underdominant threshold",
    "DESeq2 reference level",
    "Contrasts used"
  ),
  value = c(
    "0.05",
    "0.25",
    "0.5 x log2FC(GG vs AA)",
    "mid-parent +/- 0.25 log2FC",
    "log2FC(GG vs AA) +/- 0.25 log2FC",
    "+/- 0.25 log2FC around 0",
    "log2FC(hybrid vs AA) > log2FC(GG vs AA) + 0.25",
    "log2FC(hybrid vs AA) < -0.25",
    "AA",
    "GG vs AA; AG vs AA; GA vs AA"
  ),
  stringsAsFactors = FALSE
)

out_xlsx <- file.path("output", "Table_S2.dominance_thresholds.xlsx")
writexl::write_xlsx(
  list(
    parameters     = parameters,
    classification_rules = rules,
    threshold_summary    = summary_values
  ),
  path = out_xlsx
)

message("Created: ", out_xlsx)
