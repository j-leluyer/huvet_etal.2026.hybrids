suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tibble)
})

# ── Project root ────────────────────────────────────────────────────────────
args        <- commandArgs(trailingOnly = TRUE)
script_args <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", script_args[grep("--file=", script_args)])
script_dir  <- if (length(script_path) > 0) {
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

# ── Load data ────────────────────────────────────────────────────────────────
dat_file <- "output/epst_with_deseq_dominance.tsv"
if (!file.exists(dat_file)) {
  stop("Required file not found: ", dat_file,
       "\nRun script/compute_epst_quantification.R first.")
}

dat <- read.table(dat_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)

perm_file <- "output/epst_permutation_null_quantiles.tsv"
if (!file.exists(perm_file)) {
  stop("Required file not found: ", perm_file,
       "\nRun script/compute_epst_quantification.R first.")
}
perm_tbl     <- read.table(perm_file, header = TRUE, sep = "\t")
perm_thr975  <- perm_tbl$null_ePST[perm_tbl$quantile == 0.975]

# ── Derived columns ──────────────────────────────────────────────────────────
dat <- dat %>%
  filter(!is.na(ePST), ePST > 0, !is.na(log2FoldChange)) %>%
  mutate(
    neg_log10_ePST  = -log10(ePST),
    is_outlier      = ePST > perm_thr975,
    class_AG        = factor(class_AG,
                             levels = c("Ambiguous", "Additive",
                                        "AA_dominant", "GG_dominant",
                                        "Overdominant", "Underdominant"))
  )

# ── Colour palette (matching the rest of the paper) ─────────────────────────
dom_cols <- c(
  "Ambiguous"    = "grey80",
  "Additive"     = "#377EB8",
  "AA_dominant"  = "#4DAF4A",
  "GG_dominant"  = "#984EA3",
  "Overdominant" = "#E41A1C",
  "Underdominant"= "#FF7F00"
)

# ── Reference lines ──────────────────────────────────────────────────────────
vline_lfc      <- c(-2, 2)                      # user-requested vertical lines
hline_threshold <- -log10(perm_thr975)          # ePST permutation q97.5 in −log10 scale

# Split background from coloured points for rendering speed / legibility
dat_bg  <- filter(dat, class_AG == "Ambiguous")
dat_col <- filter(dat, class_AG != "Ambiguous")

# ── Volcano plot ─────────────────────────────────────────────────────────────
p_volcano <- ggplot(dat, aes(x = log2FoldChange, y = neg_log10_ePST)) +

  # background (ambiguous) layer
  geom_point(data = dat_bg,
             color = "grey82", alpha = 0.35, size = 0.6, shape = 16) +

  # coloured dominance-class layer
  geom_point(data = dat_col,
             aes(color = class_AG),
             alpha = 0.65, size = 0.75, shape = 16) +

  # highlight ePST outliers with a black ring
  geom_point(data = filter(dat, is_outlier),
             shape = 21, size = 1.4,
             fill  = NA, color = "black", stroke = 0.45) +

  # vertical lines at ±2 log2FC
  geom_vline(xintercept = vline_lfc,
             linetype = "dashed", linewidth = 0.5, color = "grey30") +

  # horizontal reference at ePST permutation threshold
  geom_hline(yintercept = hline_threshold,
             linetype   = "dotted", linewidth = 0.5, color = "firebrick") +

  scale_color_manual(values = dom_cols, drop = FALSE,
                     name   = "AG dominance class") +

  labs(
    x     = expression(log[2]*"FC  (GG vs AA)"),
    y     = expression(-log[10]*"(ePST)"),
    title = "Volcano plot: expression divergence (ePST) vs. fold-change",
    subtitle = paste0(
      "Dashed lines: |log2FC| = 2  |  ",
      "Dotted red line: permutation q97.5 threshold (ePST = ",
      round(perm_thr975, 3), ")"
    )
  ) +

  theme_classic(base_size = 11) +
  theme(
    legend.position   = "right",
    legend.key.size   = unit(0.45, "cm"),
    plot.title        = element_text(size = 11, face = "bold"),
    plot.subtitle     = element_text(size = 8, color = "grey40")
  ) +

  guides(color = guide_legend(override.aes = list(size = 2.5, alpha = 1)))

# ── Save ─────────────────────────────────────────────────────────────────────
ggsave("output/figS_volcano_epst.png", p_volcano,
       width = 7, height = 5.5, dpi = 300)
ggsave("output/figS_volcano_epst.pdf", p_volcano,
       width = 7, height = 5.5)

message("Saved:\n  output/figS_volcano_epst.png\n  output/figS_volcano_epst.pdf")
