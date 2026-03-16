#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
  library(patchwork)
})

# Uses existing hybrid network modules from:
# output/networkConstruction_hybrids_sft10.Rda
# and tests preservation AG (reference) -> GA (test).

cor <- WGCNA::cor
bicor <- WGCNA::bicor

args <- commandArgs(trailingOnly = TRUE)
script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", script_args[grep("--file=", script_args)])
script_dir <- if (length(script_path) > 0) {
  dirname(normalizePath(script_path[1], mustWork = TRUE))
} else {
  getwd()
}

project_dir <- if (length(args) >= 1) {
  normalizePath(args[1], mustWork = TRUE)
} else {
  normalizePath(file.path(script_dir, ".."), mustWork = TRUE)
}

setwd(project_dir)
dir.create("output", showWarnings = FALSE, recursive = TRUE)

load("output/dataInput_hybrids_ortho.Rda")
load("output/networkConstruction_hybrids_sft10.Rda")

if (exists("coldata") && "Cross" %in% colnames(coldata)) {
  coldata_use <- as.data.frame(coldata)
} else {
  coldata_use <- read.table("data/metadata.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  rownames(coldata_use) <- coldata_use$Ind
}

if (!"Cross" %in% colnames(coldata_use)) {
  stop("`Cross` column not found in metadata.")
}

datExpr_all <- datExpr
common_samples <- intersect(rownames(datExpr_all), rownames(coldata_use))
datExpr_all <- datExpr_all[common_samples, , drop = FALSE]
coldata_use <- coldata_use[common_samples, , drop = FALSE]

ag_samples <- rownames(coldata_use)[coldata_use$Cross == "AG"]
ga_samples <- rownames(coldata_use)[coldata_use$Cross == "GA"]

stopifnot(length(ag_samples) >= 3, length(ga_samples) >= 3)

if (exists("moduleColors_all")) {
  moduleColors_global <- moduleColors_all
} else if (exists("moduleColors")) {
  moduleColors_global <- moduleColors
} else {
  stop("No module color vector found in networkConstruction_hybrids_sft10.Rda")
}

if (is.null(names(moduleColors_global))) {
  names(moduleColors_global) <- colnames(datExpr_all)
}

# Split data and align genes
expr_ag <- datExpr_all[ag_samples, , drop = FALSE]
expr_ga <- datExpr_all[ga_samples, , drop = FALSE]

commonGenes <- intersect(colnames(expr_ag), colnames(expr_ga))
commonGenes <- commonGenes[match(commonGenes, colnames(expr_ag))]

expr_ag <- expr_ag[, commonGenes, drop = FALSE]
expr_ga <- expr_ga[, commonGenes, drop = FALSE]

moduleColors_use <- moduleColors_global[commonGenes]
keep <- !is.na(moduleColors_use) & moduleColors_use != "grey"
moduleColors_use <- moduleColors_use[keep]

expr_ag <- expr_ag[, names(moduleColors_use), drop = FALSE]
expr_ga <- expr_ga[, names(moduleColors_use), drop = FALSE]

stopifnot(identical(colnames(expr_ag), colnames(expr_ga)))

multiExpr <- list(
  REF = list(data = expr_ag),
  TEST = list(data = expr_ga)
)
colorList <- list(REF = moduleColors_use)

pres <- WGCNA::modulePreservation(
  multiExpr,
  colorList,
  referenceNetworks = 1,
  nPermutations = 1000,
  randomSeed = 123,
  verbose = 3
)

zTab <- pres$preservation$Z[[1]][[1]]
obsTab <- pres$preservation$observed[[1]][[1]]

get_Zsummary <- function(z) {
  if ("Zsummary.pres" %in% colnames(z)) return(z[, "Zsummary.pres"])
  if ("Zsummary" %in% colnames(z)) return(z[, "Zsummary"])
  dens <- intersect(colnames(z), c("Zdensity.pres", "Zdensity"))
  conn <- intersect(colnames(z), c("Zconnectivity.pres", "Zconnectivity"))
  if (length(dens) && length(conn)) return((z[, dens[1]] + z[, conn[1]]) / 2)
  stop("No Zsummary column found in preservation output.")
}

Zsummary <- get_Zsummary(zTab)
medianRank <- if ("medianRank.pres" %in% colnames(obsTab)) {
  obsTab[, "medianRank.pres"]
} else if ("medianRank" %in% colnames(obsTab)) {
  obsTab[, "medianRank"]
} else {
  NA_real_
}

summary_tab <- data.frame(
  module = rownames(zTab),
  moduleSize = zTab[, "moduleSize"],
  Zsummary = Zsummary,
  medianRank = medianRank,
  comparison = "AG_to_GA",
  stringsAsFactors = FALSE
)
summary_tab <- summary_tab[summary_tab$module != "grey", , drop = FALSE]
summary_tab <- summary_tab[order(-summary_tab$Zsummary), , drop = FALSE]

write.table(
  summary_tab,
  file = "output/module_preservation_summary_AG_to_GA.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

df <- summary_tab
df$module <- factor(df$module, levels = rev(df$module))
df$presClass <- cut(df$Zsummary, breaks = c(-Inf, 2, 10, Inf), labels = c("None", "Moderate", "Strong"))

p_z <- ggplot(df, aes(module, Zsummary)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = c(2, 10), linetype = "dashed") +
  coord_flip() +
  labs(title = "Module preservation (AG -> GA)", x = NULL, y = "Zsummary") +
  theme_classic()

p_rank <- ggplot(df, aes(module, medianRank)) +
  geom_col(fill = "grey60") +
  coord_flip() +
  labs(title = "Median rank (AG -> GA)", x = NULL, y = "medianRank") +
  theme_classic()

p_combo <- p_z + p_rank

ggsave("output/module_preservation_zsummary_AG_to_GA.png", p_z, width = 8, height = 6, dpi = 300)
ggsave("output/module_preservation_medianRank_AG_to_GA.png", p_rank, width = 8, height = 6, dpi = 300)
ggsave("output/module_preservation_combined_AG_to_GA.png", p_combo, width = 14, height = 6, dpi = 300)

message("Created hybrid preservation outputs (AG -> GA).")
