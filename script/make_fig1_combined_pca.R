#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
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
vsd.fast <- vst(dds, fitType = "local", blind = TRUE)

theme_glob <- theme(
  axis.text.x = element_text(colour = "black", size = 12),
  axis.text.y = element_text(colour = "black", size = 12),
  axis.title.x = element_text(colour = "black", size = 12, face = "bold"),
  axis.title.y = element_text(colour = "black", size = 12, face = "bold"),
  panel.background = element_blank(),
  panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.5),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  plot.margin = unit(c(1, 1, 1, 1), "line"),
  legend.title = element_blank(),
  aspect.ratio = 1,
  legend.background = element_blank(),
  legend.key = element_rect(fill = NA),
  legend.text = element_text(colour = "black", size = 12, face = "bold")
)

couleurs <- c("dodgerblue2", "tomato2", "orange2", "black")

keep_nuc <- !(rownames(vsd.fast) %in% mito_ids)
plotDataNuc <- plotPCA(vsd.fast[keep_nuc, ], intgroup = c("Cross", "gi", "weight"), returnData = TRUE)
percentVarNuc <- round(100 * attr(plotDataNuc, "percentVar"))
PCANuc <- ggplot(plotDataNuc, aes(PC1, PC2, fill = Cross)) +
  geom_point(size = 3, shape = 21, alpha = 0.7, stroke = 1) +
  xlab(paste0("PC1: ", percentVarNuc[1], "% variance")) +
  ylab(paste0("PC2: ", percentVarNuc[2], "% variance")) +
  coord_fixed()

graph.pcaNuc <- PCANuc + theme_glob + scale_color_manual(values = couleurs) +
  scale_fill_manual(values = c("AA" = "tomato2", "GG" = "blue", "AG" = "orange2", "GA" = "grey50"),
                    labels = c("AA" = "Cangulata", "GG" = "Cgigas", "AG" = "Cang-Cgig", "GA" = "Cgig-Cang")) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  annotate(geom = "text", x = -10, y = 30, label = "Nucl. orthogroups (N ~ 20k)", color = "black", size = 4)

keep_mito <- rownames(vsd.fast) %in% mito_ids
plotDataMito <- plotPCA(vsd.fast[keep_mito, ], intgroup = c("Cross", "gi", "weight"), returnData = TRUE)
percentVarMito <- round(100 * attr(plotDataMito, "percentVar"))
PCAMito <- ggplot(plotDataMito, aes(PC1, PC2, fill = Cross)) +
  geom_point(size = 3, shape = 21, alpha = 0.7, stroke = 1) +
  xlab(paste0("PC1: ", percentVarMito[1], "% variance")) +
  ylab(paste0("PC2: ", percentVarMito[2], "% variance")) +
  coord_fixed()

graph.pcaMito <- PCAMito + theme_glob + scale_color_manual(values = couleurs) +
  scale_fill_manual(values = c("AA" = "tomato2", "GG" = "blue", "AG" = "orange2", "GA" = "grey50"),
                    labels = c("AA" = "Cangulata", "GG" = "Cgigas", "AG" = "Cang-Cgig", "GA" = "Cgig-Cang")) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  annotate(geom = "text", x = 0, y = 3, label = "Mito. orthogroups (N = 37)", color = "black", size = 4)

combined_plot_fig1 <- (graph.pcaNuc | graph.pcaMito) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "right")

ggsave("output/fig1_combined_pca.png", combined_plot_fig1, width = 12, height = 6, dpi = 300)
ggsave("output/fig1_combined_pca.pdf", combined_plot_fig1, width = 12, height = 6)

message("Created output/fig1_combined_pca.png and output/fig1_combined_pca.pdf")
