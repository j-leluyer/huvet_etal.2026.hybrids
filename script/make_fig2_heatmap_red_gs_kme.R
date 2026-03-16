#!/usr/bin/env Rscript

suppressPackageStartupMessages({
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

cross <- as.data.frame(datTraits$Cross)
names(cross) <- "cross"

modNames <- substring(names(MEs), 3)
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "p"))
geneTraitSignificance <- as.data.frame(cor(datExpr, cross, use = "p"))

module <- "red"
column <- match(module, modNames)
moduleGenes <- moduleColors == module

if (is.na(column) || sum(moduleGenes) == 0) {
  stop("Module 'red' was not found in current network results.")
}

df_module_red <- data.frame(
  kME = abs(geneModuleMembership[moduleGenes, column]),
  GS = abs(geneTraitSignificance[moduleGenes, 1])
)

ct_red <- cor.test(df_module_red$kME, df_module_red$GS, method = "pearson")

p_module_red <- ggplot(df_module_red, aes(x = kME, y = GS)) +
  geom_point(color = module, alpha = 0.7, size = 1.5) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
  labs(
    x = paste("Module Membership in", module, "module"),
    y = "Gene significance for cross",
    title = "Red module: GS vs kME"
  ) +
  annotate(
    "text",
    x = 0.05,
    y = 0.95,
    hjust = 0,
    vjust = 1,
    label = paste0(
      "r = ", sprintf("%.2f", unname(ct_red$estimate)),
      "\np = ", format.pval(ct_red$p.value, digits = 2, eps = 1e-4)
    ),
    size = 4
  ) +
  theme_classic()

figure2 <- (p_heatmap | p_module_red) +
  plot_layout(widths = c(1.8, 1)) +
  plot_annotation(tag_levels = "A")

ggsave("output/fig2_heatmap_red_gs_kme.png", figure2, width = 16, height = 7, dpi = 300)
ggsave("output/fig2_heatmap_red_gs_kme.pdf", figure2, width = 16, height = 7)

message("Created output/fig2_heatmap_red_gs_kme.png and output/fig2_heatmap_red_gs_kme.pdf")
