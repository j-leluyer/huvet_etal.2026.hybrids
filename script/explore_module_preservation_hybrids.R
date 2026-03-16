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

# Remove genes with zero variance in either group before preservation
var_keep <- apply(expr_ag, 2, var, na.rm = TRUE) > 0 &
  apply(expr_ga, 2, var, na.rm = TRUE) > 0
expr_ag <- expr_ag[, var_keep, drop = FALSE]
expr_ga <- expr_ga[, var_keep, drop = FALSE]
moduleColors_use <- moduleColors_use[colnames(expr_ag)]

gsg_ag <- goodSamplesGenes(expr_ag, verbose = 0)
gsg_ga <- goodSamplesGenes(expr_ga, verbose = 0)
common_good_genes <- colnames(expr_ag)[gsg_ag$goodGenes & gsg_ga$goodGenes]
expr_ag <- expr_ag[gsg_ag$goodSamples, common_good_genes, drop = FALSE]
expr_ga <- expr_ga[gsg_ga$goodSamples, common_good_genes, drop = FALSE]
moduleColors_use <- moduleColors_use[common_good_genes]

stopifnot(ncol(expr_ag) > 0, ncol(expr_ga) > 0)
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
  nPermutations = 10,
  randomSeed = 123,
  verbose = 3
)

extract_first_table <- function(x) {
  if (is.data.frame(x)) return(x)
  if (is.matrix(x)) return(as.data.frame(x))
  if (is.list(x)) {
    for (i in seq_along(x)) {
      out <- extract_first_table(x[[i]])
      if (!is.null(out)) return(out)
    }
  }
  NULL
}

zTab <- extract_first_table(pres$preservation$Z)
obsTab <- extract_first_table(pres$preservation$observed)

if (is.null(zTab) && is.null(obsTab)) {
  stop("Could not extract preservation summary tables from WGCNA output.")
}

get_Zsummary <- function(z) {
  if (is.null(z) || nrow(z) == 0) return(numeric(0))
  zsum_col <- grep("^Zsummary|Zsummary", colnames(z), ignore.case = TRUE, value = TRUE)
  if (length(zsum_col) > 0) return(z[, zsum_col[1]])
  dens <- grep("Zdensity", colnames(z), ignore.case = TRUE, value = TRUE)
  conn <- grep("Zconnectivity", colnames(z), ignore.case = TRUE, value = TRUE)
  if (length(dens) && length(conn)) return((z[, dens[1]] + z[, conn[1]]) / 2)
  rep(NA_real_, nrow(z))
}

if (!is.null(zTab) && nrow(zTab) > 0) {
  module_names <- rownames(zTab)
  n_mod <- nrow(zTab)
} else {
  module_names <- rownames(obsTab)
  n_mod <- nrow(obsTab)
}

if (is.null(module_names)) {
  module_names <- paste0("module_", seq_len(n_mod))
}

Zsummary <- if (!is.null(zTab) && nrow(zTab) > 0) get_Zsummary(zTab) else rep(NA_real_, n_mod)
medianRank <- if (!is.null(obsTab) && "medianRank.pres" %in% colnames(obsTab)) {
  obsTab[, "medianRank.pres"]
} else if (!is.null(obsTab) && "medianRank" %in% colnames(obsTab)) {
  obsTab[, "medianRank"]
} else {
  rep(NA_real_, n_mod)
}

module_size_col <- if (!is.null(zTab) && "moduleSize" %in% colnames(zTab)) {
  "moduleSize"
} else {
  if (!is.null(zTab)) {
    grep("module.*size|size", colnames(zTab), ignore.case = TRUE, value = TRUE)[1]
  } else if (!is.null(obsTab)) {
    grep("module.*size|size", colnames(obsTab), ignore.case = TRUE, value = TRUE)[1]
  } else {
    NA_character_
  }
}

summary_tab <- data.frame(
  module = module_names,
  moduleSize = if (!is.na(module_size_col) && !is.null(zTab) && module_size_col %in% colnames(zTab)) {
    zTab[, module_size_col]
  } else if (!is.na(module_size_col) && !is.null(obsTab) && module_size_col %in% colnames(obsTab)) {
    obsTab[, module_size_col]
  } else {
    rep(NA_real_, n_mod)
  },
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

if (exists("MEs_all")) {
  me_samples <- intersect(rownames(MEs_all), rownames(coldata_use))
  MEs_tmp <- MEs_all[me_samples, , drop = FALSE]
  meta_tmp <- coldata_use[me_samples, , drop = FALSE]

  ag_me <- rownames(meta_tmp)[meta_tmp$Cross == "AG"]
  ga_me <- rownames(meta_tmp)[meta_tmp$Cross == "GA"]

  if (length(ag_me) >= 3 && length(ga_me) >= 3) {
    me_diff <- data.frame(
      module = colnames(MEs_tmp),
      mean_ref = colMeans(MEs_tmp[ag_me, , drop = FALSE]),
      mean_test = colMeans(MEs_tmp[ga_me, , drop = FALSE]),
      stringsAsFactors = FALSE
    )
    me_diff$delta <- me_diff$mean_ref - me_diff$mean_test
    me_diff$p <- sapply(colnames(MEs_tmp), function(me) {
      suppressWarnings(wilcox.test(MEs_tmp[ag_me, me], MEs_tmp[ga_me, me])$p.value)
    })
    me_diff$fdr <- p.adjust(me_diff$p, "fdr")
    me_diff <- me_diff[order(me_diff$fdr), , drop = FALSE]

    write.table(
      me_diff,
      file = "output/module_eigengene_diff_AG_to_GA.tsv",
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
}

message("Created hybrid preservation outputs (AG -> GA).")
