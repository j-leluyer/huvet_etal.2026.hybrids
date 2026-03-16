#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(WGCNA)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
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

load("output/dataInput_hybrids_ortho.Rda")
load("output/networkConstruction_hybrids_sft10.Rda")
load("output/step1.hybrids.ortho.Rda")

MEs <- orderMEs(MEs)

datTraits_plot <- datTraits[, !colnames(datTraits) %in% c("AA", "GG", "matA", "matG"), drop = FALSE]

trait_method <- ifelse(
  colnames(datTraits_plot) %in% c("gi", "weight", "lustre", "quality"),
  "pearson",
  "spearman"
)

moduleTraitCor <- matrix(
  NA,
  nrow = ncol(MEs),
  ncol = ncol(datTraits_plot),
  dimnames = list(colnames(MEs), colnames(datTraits_plot))
)

moduleTraitPvalue <- matrix(
  NA,
  nrow = ncol(MEs),
  ncol = ncol(datTraits_plot),
  dimnames = list(colnames(MEs), colnames(datTraits_plot))
)

for (i in seq_len(ncol(MEs))) {
  for (j in seq_len(ncol(datTraits_plot))) {
    x <- MEs[, i]
    y <- datTraits_plot[, j]
    ok <- complete.cases(x, y)

    if (sum(ok) >= 3) {
      tmp <- if (trait_method[j] == "spearman") {
        suppressWarnings(cor.test(x[ok], y[ok], method = "spearman", exact = FALSE))
      } else {
        suppressWarnings(cor.test(x[ok], y[ok], method = "pearson"))
      }
      moduleTraitCor[i, j] <- unname(tmp$estimate)
      moduleTraitPvalue[i, j] <- tmp$p.value
    }
  }
}
moduleTraitAdjPvalue <- matrix(
  p.adjust(as.vector(moduleTraitPvalue), method = "fdr"),
  nrow = nrow(moduleTraitPvalue),
  ncol = ncol(moduleTraitPvalue),
  dimnames = dimnames(moduleTraitPvalue)
)

cor_long <- as.data.frame(moduleTraitCor) %>%
  rownames_to_column("module") %>%
  pivot_longer(-module, names_to = "trait", values_to = "correlation")

padj_long <- as.data.frame(moduleTraitAdjPvalue) %>%
  rownames_to_column("module") %>%
  pivot_longer(-module, names_to = "trait", values_to = "adjPval")

heatmap_df <- cor_long %>%
  left_join(padj_long, by = c("module", "trait")) %>%
  mutate(
    label = ifelse(adjPval < 0.01,
                   sprintf("%.2f\n(adj.P=%.2g)", correlation, adjPval),
                   "")
  )

heatmap_df$module <- factor(heatmap_df$module, levels = rev(colnames(MEs)))
heatmap_df$trait <- factor(heatmap_df$trait, levels = colnames(datTraits_plot))
heatmap_df$module <- gsub("^ME", "", heatmap_df$module)

p_heatmap <- ggplot(heatmap_df, aes(x = trait, y = module, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 3, lineheight = 0.9) +
  scale_fill_gradient2(
    low = "#3B4CC0", mid = "white", high = "#B40426",
    midpoint = 0, limits = c(-1, 1), name = "Correlation"
  ) +
  labs(x = NULL, y = NULL, title = "Module-trait heatmap") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

dds_ag_ga <- dds[, dds$Cross %in% c("AG", "GA")]
dds_ag_ga$Cross <- droplevels(dds_ag_ga$Cross)
design(dds_ag_ga) <- ~ Cross
dds_ag_ga <- DESeq(dds_ag_ga)

res_ag_vs_ga <- results(dds_ag_ga, contrast = c("Cross", "AG", "GA"))
res_df <- as.data.frame(res_ag_vs_ga)
res_df$gene <- rownames(res_df)

write.table(
  res_df,
  file = "output/Table_S1_DESeq2_AG_vs_GA.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

res_df <- res_df[!is.na(res_df$padj) & !is.na(res_df$log2FoldChange), , drop = FALSE]

volc_all <- res_df
volc_all$contrast <- "AG vs GA"
volc_all$signif <- "NS"
volc_all$signif[volc_all$log2FoldChange >= 2 & volc_all$padj < 0.01] <- "Up"
volc_all$signif[volc_all$log2FoldChange <= -2 & volc_all$padj < 0.01] <- "Down"
volc_all$signif <- factor(volc_all$signif, levels = c("Up", "Down", "NS"))

vol_cols <- c(
  "Up" = "#B22222",
  "Down" = "#1B4F9C",
  "NS" = "grey80"
)

p_volcano <- ggplot(volc_all, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = signif), alpha = 0.6, size = 1, na.rm = TRUE) +
  scale_color_manual(values = vol_cols, drop = FALSE) +
  facet_wrap(~contrast, ncol = 1) +
  theme_bw() +
  labs(
    x = "log2 Fold Change",
    y = "-log10 adjusted p-value",
    color = "Significance"
  ) +
  geom_vline(xintercept = c(-2, 2), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.01), linetype = "dashed")

figure2 <- (p_heatmap | p_volcano) +
  plot_layout(widths = c(1.8, 1)) +
  plot_annotation(tag_levels = "A")

ggsave("output/fig2_heatmap_red_gs_kme.png", figure2, width = 16, height = 7, dpi = 300)
ggsave("output/fig2_heatmap_red_gs_kme.pdf", figure2, width = 16, height = 7)

message("Created output/fig2_heatmap_red_gs_kme.png and output/fig2_heatmap_red_gs_kme.pdf")
