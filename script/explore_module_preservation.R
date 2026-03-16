suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
  library(patchwork)
})

# Ensure WGCNA correlation functions are used
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

datExpr_all <- datExpr

if (exists("coldata") && "Cross" %in% colnames(coldata)) {
  coldata_use <- as.data.frame(coldata)
} else {
  coldata_use <- read.table(
    "data/metadata.txt",
    header = TRUE,
    sep = "\t",
    stringsAsFactors = FALSE
  )
  rownames(coldata_use) <- coldata_use$Ind
}

if (!"Cross" %in% colnames(coldata_use)) {
  stop("`Cross` column not found in metadata.")
}

common_samples <- intersect(rownames(datExpr_all), rownames(coldata_use))
if (length(common_samples) < 6) {
  stop("Too few shared samples between `datExpr` and metadata.")
}

datExpr_all <- datExpr_all[common_samples, , drop = FALSE]
coldata_use <- coldata_use[common_samples, , drop = FALSE]
coldata_use$Cross <- as.character(coldata_use$Cross)

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
    stop("Module colors are unnamed and do not match the expression gene dimension.")
  }
}

get_Zsummary <- function(z) {
  if ("Zsummary.pres" %in% colnames(z)) return(z[, "Zsummary.pres"])
  if ("Zsummary" %in% colnames(z)) return(z[, "Zsummary"])
  dens <- intersect(colnames(z), c("Zdensity.pres", "Zdensity"))
  conn <- intersect(colnames(z), c("Zconnectivity.pres", "Zconnectivity"))
  if (length(dens) && length(conn)) return((z[, dens[1]] + z[, conn[1]]) / 2)
  stop("No usable Zsummary information found.")
}

run_preservation <- function(ref_label, test_label) {
  message("\nRunning module preservation: ", ref_label, " -> ", test_label)

  ref_samples <- rownames(coldata_use)[coldata_use$Cross == ref_label]
  test_samples <- rownames(coldata_use)[coldata_use$Cross == test_label]

  if (length(ref_samples) < 3 || length(test_samples) < 3) {
    stop(
      "Insufficient samples for comparison ", ref_label, " vs ", test_label,
      " (", ref_label, "=", length(ref_samples), ", ",
      test_label, "=", length(test_samples), ")."
    )
  }

  datExpr_ref <- datExpr_all[ref_samples, , drop = FALSE]
  datExpr_test <- datExpr_all[test_samples, , drop = FALSE]

  commonGenes <- intersect(colnames(datExpr_ref), colnames(datExpr_test))
  commonGenes <- commonGenes[match(commonGenes, colnames(datExpr_ref))]

  datExpr_ref <- datExpr_ref[, commonGenes, drop = FALSE]
  datExpr_test <- datExpr_test[, commonGenes, drop = FALSE]

  moduleColors_use <- moduleColors_global[commonGenes]
  keep <- !is.na(moduleColors_use) & moduleColors_use != "grey"

  moduleColors_use <- moduleColors_use[keep]
  datExpr_ref <- datExpr_ref[, names(moduleColors_use), drop = FALSE]
  datExpr_test <- datExpr_test[, names(moduleColors_use), drop = FALSE]

  stopifnot(identical(colnames(datExpr_ref), colnames(datExpr_test)))
  stopifnot(identical(names(moduleColors_use), colnames(datExpr_ref)))

  multiExpr <- list(
    REF = list(data = datExpr_ref),
    TEST = list(data = datExpr_test)
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
    comparison = paste0(ref_label, "_to_", test_label),
    stringsAsFactors = FALSE
  )

  pres_summary <- pres_summary[pres_summary$module != "grey", , drop = FALSE]
  pres_summary <- pres_summary[order(-pres_summary$Zsummary), , drop = FALSE]

  prefix <- paste0(ref_label, "_to_", test_label)

  write.table(
    pres_summary,
    file = file.path("output", paste0("module_preservation_summary_", prefix, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  df <- pres_summary
  df$module <- factor(df$module, levels = rev(df$module))
  df$presClass <- cut(
    df$Zsummary,
    breaks = c(-Inf, 2, 10, Inf),
    labels = c("None", "Moderate", "Strong")
  )

  p_zsummary <- ggplot(df, aes(x = module, y = Zsummary)) +
    geom_col(fill = "steelblue") +
    geom_hline(yintercept = c(2, 10), linetype = "dashed") +
    coord_flip() +
    labs(
      title = paste0("Module preservation (", ref_label, " -> ", test_label, ")"),
      y = "Zsummary",
      x = NULL
    ) +
    theme_classic()

  p_class <- ggplot(df, aes(x = module, y = Zsummary, fill = presClass)) +
    geom_col() +
    geom_hline(yintercept = c(2, 10), linetype = "dashed") +
    coord_flip() +
    scale_fill_manual(values = c("None" = "grey70", "Moderate" = "orange", "Strong" = "steelblue")) +
    labs(
      title = paste0("Preservation class (", ref_label, " -> ", test_label, ")"),
      y = "Zsummary",
      x = NULL,
      fill = "Class"
    ) +
    theme_classic()

  df2 <- pres_summary[order(pres_summary$medianRank), , drop = FALSE]
  df2$module <- factor(df2$module, levels = rev(df2$module))

  p_median <- ggplot(df2, aes(x = module, y = medianRank)) +
    geom_col(fill = "darkgrey") +
    coord_flip() +
    labs(
      title = paste0("Median rank (", ref_label, " -> ", test_label, ")"),
      y = "medianRank (lower = more preserved)",
      x = NULL
    ) +
    theme_classic()

  p_combo <- p_zsummary + p_median

  ggsave(file.path("output", paste0("module_preservation_zsummary_", prefix, ".png")), p_zsummary, width = 8, height = 6, dpi = 300)
  ggsave(file.path("output", paste0("module_preservation_class_", prefix, ".png")), p_class, width = 8, height = 6, dpi = 300)
  ggsave(file.path("output", paste0("module_preservation_medianRank_", prefix, ".png")), p_median, width = 8, height = 6, dpi = 300)
  ggsave(file.path("output", paste0("module_preservation_combined_", prefix, ".png")), p_combo, width = 14, height = 6, dpi = 300)

  if (exists("MEs_all")) {
    MEs_use <- MEs_all
  } else if (exists("MEs")) {
    MEs_use <- MEs
  } else {
    MEs_use <- NULL
  }

  if (!is.null(MEs_use)) {
    ref_me <- intersect(rownames(MEs_use), ref_samples)
    test_me <- intersect(rownames(MEs_use), test_samples)

    if (length(ref_me) >= 3 && length(test_me) >= 3) {
      ME_diff <- data.frame(
        module = colnames(MEs_use),
        mean_ref = colMeans(MEs_use[ref_me, , drop = FALSE]),
        mean_test = colMeans(MEs_use[test_me, , drop = FALSE]),
        stringsAsFactors = FALSE
      )
      ME_diff$delta <- ME_diff$mean_ref - ME_diff$mean_test
      ME_diff$p <- sapply(colnames(MEs_use), function(me) {
        suppressWarnings(wilcox.test(MEs_use[ref_me, me], MEs_use[test_me, me])$p.value)
      })
      ME_diff$fdr <- p.adjust(ME_diff$p, "fdr")
      ME_diff$comparison <- prefix
      ME_diff <- ME_diff[order(ME_diff$fdr), , drop = FALSE]

      write.table(
        ME_diff,
        file = file.path("output", paste0("module_eigengene_diff_", prefix, ".tsv")),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }
  }

  pres_summary
}

# Sequential module-preservation workflow requested:
# 1) pure strains GG vs AA
# 2) hybrids AG vs GA
res_pure <- run_preservation("GG", "AA")
res_hybrid <- run_preservation("AG", "GA")

all_summary <- rbind(res_pure, res_hybrid)
write.table(
  all_summary,
  file = "output/module_preservation_summary_all_comparisons.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Created sequential module-preservation outputs for GG->AA and AG->GA.")
