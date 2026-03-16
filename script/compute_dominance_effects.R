suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(pheatmap)
  library(patchwork)
})

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

load("output/step1.hybrids.ortho.Rda")

if (!exists("dds") || !exists("mito_ids")) {
  stop("Missing required objects in output/step1.hybrids.ortho.Rda (need `dds` and `mito_ids`).")
}

classify_dominance <- function(res_hyb_vs_AA, res_GG_vs_AA, alpha = 0.05, delta = 0.25) {
  cls <- rep("Ambiguous", nrow(res_hyb_vs_AA))

  # Additive
  cls[
    res_hyb_vs_AA$padj < alpha &
      res_GG_vs_AA$padj < alpha &
      abs(res_hyb_vs_AA$log2FoldChange - 0.5 * res_GG_vs_AA$log2FoldChange) < delta
  ] <- "Additive"

  # GG-dominant
  cls[
    res_hyb_vs_AA$padj < alpha &
      abs(res_hyb_vs_AA$log2FoldChange - res_GG_vs_AA$log2FoldChange) < delta
  ] <- "GG_dominant"

  # AA-dominant
  cls[
    res_hyb_vs_AA$padj < alpha &
      abs(res_hyb_vs_AA$log2FoldChange) < delta
  ] <- "AA_dominant"

  # Overdominant
  cls[
    res_hyb_vs_AA$padj < alpha &
      res_GG_vs_AA$padj < alpha &
      res_hyb_vs_AA$log2FoldChange > res_GG_vs_AA$log2FoldChange + delta
  ] <- "Overdominant"

  # Underdominant
  cls[
    res_hyb_vs_AA$padj < alpha &
      res_GG_vs_AA$padj < alpha &
      res_hyb_vs_AA$log2FoldChange < -delta
  ] <- "Underdominant"

  cls
}

run_dominance <- function(dds_sub, label, alpha = 0.05, delta = 0.25) {
  design(dds_sub) <- ~ Cross
  dds_sub <- DESeq(dds_sub)

  res_GG_AA <- results(dds_sub, contrast = c("Cross", "GG", "AA"))
  res_AG_AA <- results(dds_sub, contrast = c("Cross", "AG", "AA"))
  res_GA_AA <- results(dds_sub, contrast = c("Cross", "GA", "AA"))

  denom <- 0.5 * res_GG_AA$log2FoldChange
  denom_safe <- ifelse(abs(denom) < 1e-12, NA_real_, denom)

  D_AG <- (res_AG_AA$log2FoldChange - 0.5 * res_GG_AA$log2FoldChange) / denom_safe
  D_GA <- (res_GA_AA$log2FoldChange - 0.5 * res_GG_AA$log2FoldChange) / denom_safe

  class_AG <- classify_dominance(res_AG_AA, res_GG_AA, alpha = alpha, delta = delta)
  class_GA <- classify_dominance(res_GA_AA, res_GG_AA, alpha = alpha, delta = delta)

  tab_AG <- table(class_AG)
  tab_GA <- table(class_GA)
  tab_cross <- table(class_AG, class_GA)

  asymmetric <- class_AG != class_GA & !(class_AG == "Ambiguous" & class_GA == "Ambiguous")

  res_tbl <- tibble(
    gene = rownames(res_GG_AA),
    log2FC_GG_vs_AA = res_GG_AA$log2FoldChange,
    padj_GG_vs_AA = res_GG_AA$padj,
    log2FC_AG_vs_AA = res_AG_AA$log2FoldChange,
    padj_AG_vs_AA = res_AG_AA$padj,
    log2FC_GA_vs_AA = res_GA_AA$log2FoldChange,
    padj_GA_vs_AA = res_GA_AA$padj,
    D_AG = D_AG,
    D_GA = D_GA,
    class_AG = class_AG,
    class_GA = class_GA,
    asymmetric = asymmetric
  )

  write.table(
    res_tbl,
    file = file.path("output", paste0("dominance_gene_level_", label, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  summary_tbl <- tibble(
    class = names(tab_AG),
    AG_n = as.integer(tab_AG),
    AG_prop = as.numeric(prop.table(tab_AG)[names(tab_AG)]),
    GA_n = as.integer(tab_GA[names(tab_AG)]),
    GA_prop = as.numeric(prop.table(tab_GA)[names(tab_AG)]),
    dataset = label
  )

  write.table(
    summary_tbl,
    file = file.path("output", paste0("dominance_class_summary_", label, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  cross_df <- as.data.frame(tab_cross)
  names(cross_df) <- c("class_AG", "class_GA", "n")
  write.table(
    cross_df,
    file = file.path("output", paste0("dominance_cross_table_", label, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  asym_df <- tibble(
    dataset = label,
    n_asymmetric = sum(asymmetric, na.rm = TRUE),
    frac_asymmetric = mean(asymmetric, na.rm = TRUE)
  )
  write.table(
    asym_df,
    file = file.path("output", paste0("dominance_asymmetry_summary_", label, ".tsv")),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  # Heatmap of AG vs GA class transitions
  mat <- as.matrix(tab_cross)
  pheatmap(
    mat / pmax(rowSums(mat), 1),
    display_numbers = mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    filename = file.path("output", paste0("dominance_class_heatmap_", label, ".png")),
    width = 7,
    height = 6
  )

  heat_df <- cross_df %>%
    group_by(class_AG) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()

  p_heat <- ggplot(heat_df, aes(x = class_GA, y = class_AG, fill = prop)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = n), size = 3) +
    scale_fill_gradient(low = "white", high = "steelblue", name = "Row proportion") +
    theme_classic() +
    labs(
      title = paste0("Class transitions (", label, ")"),
      x = "GA class",
      y = "AG class"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Dominance scatter for AG
  scatter_df <- tibble(
    x = res_GG_AA$log2FoldChange,
    y = res_AG_AA$log2FoldChange,
    class = class_AG
  )

  p_scatter <- ggplot(scatter_df, aes(x = x, y = y, color = class)) +
    geom_point(alpha = 0.5, size = 0.8, na.rm = TRUE) +
    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
    geom_abline(intercept = 0, slope = 0.5, color = "blue", linetype = "dashed") +
    geom_abline(intercept = 0, slope = 0, color = "darkgreen", linetype = "dashed") +
    theme_classic() +
    labs(
      title = paste0("Dominance map (", label, ", AG vs AA vs GG vs AA)"),
      x = "GG vs AA (log2FC)",
      y = "AG vs AA (log2FC)",
      color = "Class"
    )

  ggsave(
    filename = file.path("output", paste0("dominance_scatter_AG_", label, ".png")),
    plot = p_scatter,
    width = 7,
    height = 6,
    dpi = 300
  )

  message("[", label, "] Class AG:")
  print(tab_AG)
  message("[", label, "] Class GA:")
  print(tab_GA)
  message("[", label, "] AG vs GA table:")
  print(tab_cross)
  message("[", label, "] Asymmetric fraction: ", round(mean(asymmetric, na.rm = TRUE), 4))

  invisible(list(
    dds = dds_sub,
    res_GG_AA = res_GG_AA,
    res_AG_AA = res_AG_AA,
    res_GA_AA = res_GA_AA,
    D_AG = D_AG,
    D_GA = D_GA,
    class_AG = class_AG,
    class_GA = class_GA,
    asym = asymmetric,
    p_heat = p_heat,
    p_scatter = p_scatter
  ))
}

# Nuclear genes
dds_hyb <- dds[!rownames(dds) %in% mito_ids, ]
res_nuclear <- run_dominance(dds_hyb, label = "nuclear")

# Mitochondrial genes
dds_mit <- dds[rownames(dds) %in% mito_ids, ]
res_mito <- run_dominance(dds_mit, label = "mito")

# Supporting Figure 3: dominance overview
fig3 <-
  (res_nuclear$p_heat + res_nuclear$p_scatter) /
  (res_mito$p_heat + res_mito$p_scatter) +
  plot_annotation(tag_levels = "A", title = "Figure 3. Dominance effects overview")

ggsave(
  filename = "output/fig3_dominance_supporting.png",
  plot = fig3,
  width = 14,
  height = 10,
  dpi = 300
)

ggsave(
  filename = "output/fig3_dominance_supporting.pdf",
  plot = fig3,
  width = 14,
  height = 10
)

# Optional gene-level expression plots if data objects are present
if (exists("vst.mat") && exists("coldata")) {
  data.gene.merged <- as.data.frame(merge(t(vst.mat), coldata, by = 0))
  genes_to_plot <- c("ORTHO_015690", "gene-LOC128169801")

  for (g in genes_to_plot) {
    if (g %in% colnames(data.gene.merged)) {
      p <- ggplot(data.gene.merged, aes(x = Cross, y = .data[[g]], fill = as.factor(Cross))) +
        geom_point(size = 3, shape = 21, alpha = 0.6, position = position_jitter(width = 0.15, height = 0)) +
        theme_classic() +
        labs(title = g, x = "Cross", y = "VST expression", fill = "Cross")

      ggsave(
        filename = file.path("output", paste0("dominance_gene_plot_", gsub("[^A-Za-z0-9_]+", "_", g), ".png")),
        plot = p,
        width = 6,
        height = 4,
        dpi = 300
      )
    }
  }
}

message("Created dominance outputs in output/ for nuclear and mito datasets.")
