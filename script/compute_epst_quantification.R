suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tibble)
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

# Recompute dominance classes for AG (used in ePST overlay plots)
classify_dominance <- function(res_hyb_vs_AA, res_GG_vs_AA, alpha = 0.05, delta = 0.25) {
  cls <- rep("Ambiguous", nrow(res_hyb_vs_AA))

  cls[
    res_hyb_vs_AA$padj < alpha &
      res_GG_vs_AA$padj < alpha &
      abs(res_hyb_vs_AA$log2FoldChange - 0.5 * res_GG_vs_AA$log2FoldChange) < delta
  ] <- "Additive"

  cls[
    res_hyb_vs_AA$padj < alpha &
      abs(res_hyb_vs_AA$log2FoldChange - res_GG_vs_AA$log2FoldChange) < delta
  ] <- "GG_dominant"

  cls[
    res_hyb_vs_AA$padj < alpha &
      abs(res_hyb_vs_AA$log2FoldChange) < delta
  ] <- "AA_dominant"

  cls[
    res_hyb_vs_AA$padj < alpha &
      res_GG_vs_AA$padj < alpha &
      res_hyb_vs_AA$log2FoldChange > res_GG_vs_AA$log2FoldChange + delta
  ] <- "Overdominant"

  cls[
    res_hyb_vs_AA$padj < alpha &
      res_GG_vs_AA$padj < alpha &
      res_hyb_vs_AA$log2FoldChange < -delta
  ] <- "Underdominant"

  cls
}

# Nuclear only (as in initial workflow)
dds_hyb <- dds[!rownames(dds) %in% mito_ids, ]
design(dds_hyb) <- ~ Cross
dds_hyb <- DESeq(dds_hyb)

res_GG_AA <- results(dds_hyb, contrast = c("Cross", "GG", "AA"))
res_AG_AA <- results(dds_hyb, contrast = c("Cross", "AG", "AA"))
class_AG <- classify_dominance(res_AG_AA, res_GG_AA)
names(class_AG) <- rownames(res_AG_AA)

# ePST computation (AA vs GG)
vsd.fast <- vst(dds_hyb, fitType = "local", blind = TRUE)
meta <- as.data.frame(colData(vsd.fast))
expr <- assay(vsd.fast)

keep <- meta$Cross %in% c("AA", "GG")
expr2 <- expr[, keep, drop = FALSE]
meta2 <- meta[keep, , drop = FALSE]
group <- factor(meta2$Cross)

compute_ePST <- function(x, group) {
  means <- tapply(x, group, mean)
  vars <- tapply(x, group, var)

  V_between <- var(means)
  V_within <- mean(vars, na.rm = TRUE)

  V_between / (V_between + 2 * V_within)
}

ePST <- apply(expr2, 1, compute_ePST, group = group)

thr95 <- quantile(ePST, 0.95, na.rm = TRUE)
thr99 <- quantile(ePST, 0.99, na.rm = TRUE)

stats_tbl <- tibble(
  n_genes = length(ePST),
  mean_ePST = mean(ePST, na.rm = TRUE),
  sd_ePST = sd(ePST, na.rm = TRUE),
  median_ePST = quantile(ePST, 0.5, na.rm = TRUE),
  q90_ePST = quantile(ePST, 0.9, na.rm = TRUE),
  q95_ePST = thr95,
  q99_ePST = thr99,
  n_above_q95 = sum(ePST > thr95, na.rm = TRUE),
  prop_above_0_12 = mean(ePST > 0.12, na.rm = TRUE)
)

write.table(
  stats_tbl,
  file = "output/epst_summary_stats.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

epst_tbl <- tibble(
  gene = names(ePST),
  ePST = as.numeric(ePST)
) %>%
  mutate(
    epst_class = case_when(
      ePST < 0.05 ~ "Mostly plastic / shared regulation",
      ePST < 0.12 ~ "Weak differentiation",
      ePST < 0.30 ~ "Regulatory divergence",
      TRUE ~ "Strong divergence"
    )
  )

write.table(
  epst_tbl,
  file = "output/epst_gene_values.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

res_tbl <- as.data.frame(res_GG_AA)
res_tbl$gene <- rownames(res_tbl)
res_tbl$ePST <- ePST[res_tbl$gene]
res_tbl$class_AG <- class_AG[res_tbl$gene]
res_tbl$absLFC <- abs(res_tbl$log2FoldChange)
res_tbl <- res_tbl[complete.cases(res_tbl$ePST, res_tbl$log2FoldChange), , drop = FALSE]

set_strong <- with(res_tbl, padj < 1e-4 & absLFC > 2 & ePST > thr95)
set_consistent <- with(res_tbl, ePST > thr95 & absLFC < 1)
set_noisy <- with(res_tbl, absLFC > 2 & ePST < quantile(ePST, 0.5, na.rm = TRUE))

res_tbl$set_strong <- set_strong
res_tbl$set_consistent <- set_consistent
res_tbl$set_noisy <- set_noisy
res_tbl$hit <- with(res_tbl, ePST > thr95 & absLFC > 2 & padj < 1e-4)

write.table(
  res_tbl,
  file = "output/epst_with_deseq_dominance.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  subset(res_tbl, set_strong),
  file = "output/epst_set_strong.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  subset(res_tbl, set_consistent),
  file = "output/epst_set_consistent.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

write.table(
  subset(res_tbl, set_noisy),
  file = "output/epst_set_noisy.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

class_levels <- c("Ambiguous", "Additive", "AA_dominant", "GG_dominant", "Overdominant", "Underdominant")
res_tbl$class_AG <- factor(res_tbl$class_AG, levels = class_levels)

dom_cols <- c(
  "Additive" = "#377EB8",
  "AA_dominant" = "#4DAF4A",
  "GG_dominant" = "#984EA3",
  "Overdominant" = "#E41A1C",
  "Underdominant" = "#FF7F00"
)

p_box <- ggplot(res_tbl, aes(x = class_AG, y = ePST)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90") +
  geom_jitter(width = 0.2, alpha = 0.15, size = 0.5, color = "grey45") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(x = "AG dominance class", y = "ePST", title = "ePST by dominance class")

p_lfc_epst <- ggplot(res_tbl, aes(x = absLFC, y = ePST)) +
  geom_point(data = subset(res_tbl, class_AG == "Ambiguous"), color = "grey85", alpha = 0.3, size = 0.7) +
  geom_point(
    data = subset(res_tbl, class_AG != "Ambiguous"),
    aes(color = class_AG),
    alpha = 0.6,
    size = 0.8
  ) +
  scale_color_manual(values = dom_cols, drop = FALSE) +
  geom_hline(yintercept = thr95, linetype = "dashed") +
  geom_vline(xintercept = 2, linetype = "dashed") +
  theme_classic() +
  labs(x = "|log2FC| (GG vs AA)", y = "ePST", color = "AG class", title = "ePST vs expression shift")

p_hit <- ggplot(res_tbl, aes(x = absLFC, y = ePST)) +
  geom_point(aes(color = class_AG), size = 0.5, alpha = 0.25) +
  geom_point(data = subset(res_tbl, hit), size = 1.1, color = "black") +
  geom_hline(yintercept = thr95, linetype = "dashed") +
  geom_vline(xintercept = 2, linetype = "dashed") +
  theme_classic() +
  labs(x = "|log2FC| (GG vs AA)", y = "ePST", color = "AG class", title = "Top candidates highlighted")

p_bin2d <- ggplot(res_tbl, aes(x = absLFC, y = ePST)) +
  geom_bin2d(bins = 80) +
  facet_wrap(~class_AG) +
  theme_classic() +
  labs(x = "|log2FC| (GG vs AA)", y = "ePST", title = "Density by dominance class")

fig4 <- (p_box | p_lfc_epst) / (p_hit | p_bin2d) +
  plot_annotation(tag_levels = "A", title = "Figure 4. ePST quantification and dominance overlay")

ggsave("output/fig4_epst_supporting.png", fig4, width = 14, height = 10, dpi = 300)
ggsave("output/fig4_epst_supporting.pdf", fig4, width = 14, height = 10)

ggsave("output/epst_boxplot_by_class.png", p_box, width = 7, height = 5, dpi = 300)
ggsave("output/epst_vs_abslfc_colored.png", p_lfc_epst, width = 7, height = 5, dpi = 300)
ggsave("output/epst_hits_overlay.png", p_hit, width = 7, height = 5, dpi = 300)
ggsave("output/epst_bin2d_by_class.png", p_bin2d, width = 9, height = 6, dpi = 300)

message("Created ePST outputs in output/:\n",
        " - epst_summary_stats.tsv\n",
        " - epst_gene_values.tsv\n",
        " - epst_with_deseq_dominance.tsv\n",
        " - epst_set_strong.tsv\n",
        " - epst_set_consistent.tsv\n",
        " - epst_set_noisy.tsv\n",
        " - fig4_epst_supporting.png/.pdf and component figures")
