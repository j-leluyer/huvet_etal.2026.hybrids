#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(WGCNA)
  library(ggplot2)
  library(patchwork)
})

# Builds a NEW pure-strain network (AA+GG only), then tests preservation GG -> AA.
# Network settings mirror the hybrid construction basis (signed, softPower=10,
# minModuleSize=50, mergeCutHeight=0.25, Pearson correlation).

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

if (!file.exists("output/step1.hybrids.ortho.Rda")) {
  stop("Missing output/step1.hybrids.ortho.Rda. Run script/rnaseq_ortho.R first.")
}
load("output/step1.hybrids.ortho.Rda")

if (!exists("dds") || !exists("mito_ids")) {
  stop("`dds` and/or `mito_ids` not found in output/step1.hybrids.ortho.Rda")
}

# 1) Build pure-only expression matrix (samples x genes)
pure_samples <- grep("^(AA|GG)", colnames(dds), value = TRUE)
if (length(pure_samples) < 6) {
  stop("Too few AA/GG samples available for pure network construction.")
}

dds_pure <- dds[!rownames(dds) %in% mito_ids, pure_samples]
dds_pure <- estimateSizeFactors(dds_pure)
vsd_pure <- vst(dds_pure, fitType = "local", blind = TRUE)

expr0 <- as.data.frame(t(assay(vsd_pure)))

# Basic QC as in the hybrid workflow
gsg <- goodSamplesGenes(expr0, verbose = 3)
if (!gsg$allOK) {
  expr0 <- expr0[gsg$goodSamples, gsg$goodGenes, drop = FALSE]
}

# Variance filter (same threshold basis)
vari <- apply(as.matrix(expr0), 2, var, na.rm = TRUE)
keep_var <- vari > 0.05
expr_pure <- expr0[, keep_var, drop = FALSE]

if (ncol(expr_pure) < 1000) {
  message("Warning: low number of genes after filtering (", ncol(expr_pure), ").")
}

# 2) Build a NEW pure network with same settings as hybrids
softPower <- 10

net_pure <- blockwiseModules(
  expr_pure,
  power = softPower,
  TOMType = "signed",
  networkType = "signed",
  minModuleSize = 50,
  mergeCutHeight = 0.25,
  numericLabels = FALSE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3,
  corType = "pearson",
  corFnc = WGCNA::cor,
  corOptions = list(use = "p")
)

moduleColors_pure <- net_pure$colors
names(moduleColors_pure) <- colnames(expr_pure)
MEs_pure <- orderMEs(net_pure$MEs)

save(
  expr_pure,
  moduleColors_pure,
  MEs_pure,
  net_pure,
  file = "output/networkConstruction_pure_sft10.Rda"
)

# 3) Preservation GG -> AA using this NEW pure network
coldata_use <- as.data.frame(colData(dds_pure))
coldata_use$Cross <- as.character(coldata_use$Cross)
rownames(coldata_use) <- colnames(dds_pure)

common_samples <- intersect(rownames(expr_pure), rownames(coldata_use))
expr_pure <- expr_pure[common_samples, , drop = FALSE]
coldata_use <- coldata_use[common_samples, , drop = FALSE]

gg_samples <- rownames(coldata_use)[coldata_use$Cross == "GG"]
aa_samples <- rownames(coldata_use)[coldata_use$Cross == "AA"]

stopifnot(length(gg_samples) >= 3, length(aa_samples) >= 3)

expr_gg <- expr_pure[gg_samples, , drop = FALSE]
expr_aa <- expr_pure[aa_samples, , drop = FALSE]

common_genes <- intersect(colnames(expr_gg), colnames(expr_aa))
common_genes <- common_genes[match(common_genes, colnames(expr_gg))]

expr_gg <- expr_gg[, common_genes, drop = FALSE]
expr_aa <- expr_aa[, common_genes, drop = FALSE]

moduleColors_use <- moduleColors_pure[common_genes]
keep <- !is.na(moduleColors_use) & moduleColors_use != "grey"
moduleColors_use <- moduleColors_use[keep]

expr_gg <- expr_gg[, names(moduleColors_use), drop = FALSE]
expr_aa <- expr_aa[, names(moduleColors_use), drop = FALSE]

stopifnot(identical(colnames(expr_gg), colnames(expr_aa)))

multiExpr <- list(
  REF = list(data = expr_gg),
  TEST = list(data = expr_aa)
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
  comparison = "GG_to_AA",
  stringsAsFactors = FALSE
)
summary_tab <- summary_tab[summary_tab$module != "grey", , drop = FALSE]
summary_tab <- summary_tab[order(-summary_tab$Zsummary), , drop = FALSE]

write.table(
  summary_tab,
  file = "output/module_preservation_summary_GG_to_AA.tsv",
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
  labs(title = "Module preservation (GG -> AA)", x = NULL, y = "Zsummary") +
  theme_classic()

p_rank <- ggplot(df, aes(module, medianRank)) +
  geom_col(fill = "grey60") +
  coord_flip() +
  labs(title = "Median rank (GG -> AA)", x = NULL, y = "medianRank") +
  theme_classic()

p_combo <- p_z + p_rank

ggsave("output/module_preservation_zsummary_GG_to_AA.png", p_z, width = 8, height = 6, dpi = 300)
ggsave("output/module_preservation_medianRank_GG_to_AA.png", p_rank, width = 8, height = 6, dpi = 300)
ggsave("output/module_preservation_combined_GG_to_AA.png", p_combo, width = 14, height = 6, dpi = 300)

# Optional eigengene differential table (GG vs AA)
me_samples <- intersect(rownames(MEs_pure), rownames(coldata_use))
MEs_tmp <- MEs_pure[me_samples, , drop = FALSE]
meta_tmp <- coldata_use[me_samples, , drop = FALSE]

gg_me <- rownames(meta_tmp)[meta_tmp$Cross == "GG"]
aa_me <- rownames(meta_tmp)[meta_tmp$Cross == "AA"]

if (length(gg_me) >= 3 && length(aa_me) >= 3) {
  me_diff <- data.frame(
    module = colnames(MEs_tmp),
    mean_ref = colMeans(MEs_tmp[gg_me, , drop = FALSE]),
    mean_test = colMeans(MEs_tmp[aa_me, , drop = FALSE]),
    stringsAsFactors = FALSE
  )
  me_diff$delta <- me_diff$mean_ref - me_diff$mean_test
  me_diff$p <- sapply(colnames(MEs_tmp), function(me) {
    suppressWarnings(wilcox.test(MEs_tmp[gg_me, me], MEs_tmp[aa_me, me])$p.value)
  })
  me_diff$fdr <- p.adjust(me_diff$p, "fdr")
  me_diff <- me_diff[order(me_diff$fdr), , drop = FALSE]

  write.table(
    me_diff,
    file = "output/module_eigengene_diff_GG_to_AA.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
}

message("Created pure-strain network and preservation outputs (GG -> AA).")
