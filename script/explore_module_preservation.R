suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
  library(patchwork)
})

# STEP 4: MODULE PRESERVATION (WGCNA) ----
# Reference = phenotype "L" (if available), Test = phenotype "C" (if available)

# Ensure WGCNA cor is used (important if Bioconductor/NetRep loaded)
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
if (file.exists("output/step1.hybrids.ortho.Rda")) {
  load("output/step1.hybrids.ortho.Rda")
}

# datExpr from dataInput is the global expression matrix used for network build (samples x genes)
datExpr_all <- datExpr

# Build phenotype vector
if (exists("coldata") && !is.null(rownames(coldata))) {
  coldata_use <- as.data.frame(coldata)
} else {
  coldata_use <- data.frame(row.names = rownames(datExpr_all))
}

if (!"phenotype" %in% colnames(coldata_use)) {
  if ("Cross" %in% colnames(coldata_use)) {
    coldata_use$phenotype <- as.character(coldata_use$Cross)
  } else {
    stop("No `phenotype` or `Cross` column available to define groups.")
  }
}

# Keep shared samples only
common_samples <- intersect(rownames(datExpr_all), rownames(coldata_use))
if (length(common_samples) < 6) {
  stop("Too few shared samples between expression matrix and metadata.")
}
datExpr_all <- datExpr_all[common_samples, , drop = FALSE]
coldata_use <- coldata_use[common_samples, , drop = FALSE]

# Prefer L/C if present, otherwise fallback to first two phenotype levels
phen <- as.character(coldata_use$phenotype)
if (all(c("L", "C") %in% phen)) {
  ref_label <- "L"
  test_label <- "C"
} else {
  lev <- sort(unique(phen))
  if (length(lev) < 2) {
    stop("Need at least two phenotype groups for module preservation.")
  }
  ref_label <- lev[1]
  test_label <- lev[2]
  message("Using fallback phenotype groups: ", ref_label, " (reference) and ", test_label, " (test)")
}

# 4.1 Split samples
L_samples <- rownames(coldata_use)[coldata_use$phenotype == ref_label]
C_samples <- rownames(coldata_use)[coldata_use$phenotype == test_label]

if (length(L_samples) <= 10 || length(C_samples) <= 10) {
  message("Warning: one group has <= 10 samples (", ref_label, "=", length(L_samples),
          ", ", test_label, "=", length(C_samples), ").")
}
stopifnot(length(L_samples) >= 3, length(C_samples) >= 3)

# 4.2 Split expression (samples x genes)
datExpr_L <- datExpr_all[L_samples, , drop = FALSE]
datExpr_C <- datExpr_all[C_samples, , drop = FALSE]

# 4.3 Attach gene names to GLOBAL module colors and align to expression columns
if (exists("moduleColors_all")) {
  moduleColors_global <- moduleColors_all
} else if (exists("moduleColors")) {
  moduleColors_global <- moduleColors
} else {
  stop("No module color vector found (`moduleColors_all` or `moduleColors`).")
}

if (is.null(names(moduleColors_global))) {
  if (length(moduleColors_global) == ncol(datExpr_all)) {
    names(moduleColors_global) <- colnames(datExpr_all)
  } else {
    stop("Module colors are unnamed and length does not match expression columns.")
  }
}

# Genes present in both
commonGenes <- intersect(colnames(datExpr_L), colnames(datExpr_C))

# Force the SAME order in both matrices (use L's order as reference)
commonGenes <- commonGenes[match(commonGenes, colnames(datExpr_L))]

datExpr_L <- datExpr_L[, commonGenes, drop = FALSE]
datExpr_C <- datExpr_C[, commonGenes, drop = FALSE]

# subset colors by those gene names
moduleColors_global <- moduleColors_global[commonGenes]

# Optional: drop grey (recommended)
keep <- !is.na(moduleColors_global) & moduleColors_global != "grey"
moduleColors_global <- moduleColors_global[keep]
datExpr_L <- datExpr_L[, names(moduleColors_global), drop = FALSE]
datExpr_C <- datExpr_C[, names(moduleColors_global), drop = FALSE]

# Final sanity checks
stopifnot(identical(colnames(datExpr_L), colnames(datExpr_C)))
stopifnot(identical(names(moduleColors_global), colnames(datExpr_L)))

# 4.4 Build multiExpr and colorList
multiExpr <- list(
  L = list(data = datExpr_L),
  C = list(data = datExpr_C)
)

# IMPORTANT: colorList name must match the REFERENCE dataset name in multiExpr.
colorList <- list(L = moduleColors_global)

# 4.5 Run preservation
nPerm <- 1000
pres_globl <- WGCNA::modulePreservation(
  multiExpr,
  colorList,
  referenceNetworks = 1,
  nPermutations = nPerm,
  randomSeed = 123,
  verbose = 3
)

# 4.6 Extract preservation stats: reference = L, test = C
zTab <- pres_globl$preservation$Z$ref.L$inColumnsAlsoPresentIn.C
obsTab <- pres_globl$preservation$observed$ref.L$inColumnsAlsoPresentIn.C

get_Zsummary <- function(z) {
  if ("Zsummary.pres" %in% colnames(z)) return(z[, "Zsummary.pres"])
  if ("Zsummary" %in% colnames(z)) return(z[, "Zsummary"])
  dens <- intersect(colnames(z), c("Zdensity.pres", "Zdensity"))
  conn <- intersect(colnames(z), c("Zconnectivity.pres", "Zconnectivity"))
  if (length(dens) && length(conn)) return((z[, dens[1]] + z[, conn[1]]) / 2)
  stop("No Zsummary column found. Inspect colnames(zTab) to choose the right columns.")
}

Zsummary <- get_Zsummary(zTab)

medianRank <- if ("medianRank.pres" %in% colnames(obsTab)) {
  obsTab[, "medianRank.pres"]
} else if ("medianRank" %in% colnames(obsTab)) {
  obsTab[, "medianRank"]
} else {
  NA_real_
}

pres_summary <- data.frame(
  module = rownames(zTab),
  moduleSize = zTab[, "moduleSize"],
  Zsummary = Zsummary,
  medianRank = medianRank,
  stringsAsFactors = FALSE
)

pres_summary <- pres_summary[pres_summary$module != "grey", ]
pres_summary <- pres_summary[order(-pres_summary$Zsummary), ]
print(pres_summary)

write.table(
  pres_summary,
  file = "output/module_preservation_summary_L_to_C.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

df <- pres_summary
df <- df[order(df$Zsummary, decreasing = TRUE), ]
df$module <- factor(df$module, levels = df$module)

p_zsummary <- ggplot(df, aes(x = module, y = Zsummary)) +
  geom_col(fill = "steelblue") +
  geom_hline(yintercept = 2, linetype = "dashed", color = "orange") +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  coord_flip() +
  labs(title = "Module preservation (L -> C)", y = "Zsummary", x = NULL) +
  theme_classic()

df$presClass <- cut(
  df$Zsummary,
  breaks = c(-Inf, 2, 10, Inf),
  labels = c("None", "Moderate", "Strong")
)

p_class <- ggplot(df, aes(x = module, y = Zsummary, fill = presClass)) +
  geom_col() +
  geom_hline(yintercept = c(2, 10), linetype = "dashed") +
  coord_flip() +
  scale_fill_manual(values = c("None" = "grey70", "Moderate" = "orange", "Strong" = "steelblue")) +
  labs(title = "Module preservation (L reference, C test)", y = "Zsummary", x = NULL) +
  theme_classic()

df2 <- pres_summary
df2 <- df2[order(df2$medianRank), ]
df2$module <- factor(df2$module, levels = df2$module)

p_median <- ggplot(df2, aes(x = module, y = medianRank)) +
  geom_col(fill = "darkgrey") +
  coord_flip() +
  labs(title = "Module preservation (medianRank)", y = "medianRank (lower = more preserved)", x = NULL) +
  theme_classic()

p_combo <- p_zsummary + p_median

ggsave("output/module_preservation_zsummary_L_to_C.png", p_zsummary, width = 8, height = 6, dpi = 300)
ggsave("output/module_preservation_class_L_to_C.png", p_class, width = 8, height = 6, dpi = 300)
ggsave("output/module_preservation_medianRank_L_to_C.png", p_median, width = 8, height = 6, dpi = 300)
ggsave("output/module_preservation_combined_L_to_C.png", p_combo, width = 14, height = 6, dpi = 300)

# Optional: quick look at available columns
print(colnames(zTab))
print(colnames(obsTab))

# Check module eigengene properties between L and C groups
if (exists("MEs_all")) {
  MEs_use <- MEs_all
} else if (exists("MEs")) {
  MEs_use <- MEs
} else {
  stop("No module eigengenes found (`MEs_all` or `MEs`).")
}

L <- rownames(coldata_use)[coldata_use$phenotype == ref_label]
C <- rownames(coldata_use)[coldata_use$phenotype == test_label]
common_me_samples <- intersect(rownames(MEs_use), intersect(L, C))

# Keep all group-specific samples that exist in MEs
L_me <- intersect(rownames(MEs_use), L)
C_me <- intersect(rownames(MEs_use), C)
stopifnot(length(L_me) >= 3, length(C_me) >= 3)

ME_diff <- data.frame(
  module = colnames(MEs_use),
  mean_L = colMeans(MEs_use[L_me, , drop = FALSE]),
  mean_C = colMeans(MEs_use[C_me, , drop = FALSE]),
  stringsAsFactors = FALSE
)
ME_diff$delta <- ME_diff$mean_L - ME_diff$mean_C

ME_diff$p <- sapply(colnames(MEs_use), function(me) {
  suppressWarnings(wilcox.test(MEs_use[L_me, me], MEs_use[C_me, me])$p.value)
})
ME_diff$fdr <- p.adjust(ME_diff$p, "fdr")

ME_diff <- ME_diff[order(ME_diff$fdr), ]
print(head(ME_diff, 15))

write.table(
  ME_diff,
  file = "output/module_eigengene_diff_L_vs_C.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Created:\n",
        " - output/module_preservation_summary_L_to_C.tsv\n",
        " - output/module_preservation_zsummary_L_to_C.png\n",
        " - output/module_preservation_class_L_to_C.png\n",
        " - output/module_preservation_medianRank_L_to_C.png\n",
        " - output/module_preservation_combined_L_to_C.png\n",
        " - output/module_eigengene_diff_L_vs_C.tsv")
