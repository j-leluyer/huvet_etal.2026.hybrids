#!/usr/bin/Rscript

# Script hybrids - ortholog

rm(list = ls())


# Load libraries ----
required_packages <- c(
  "vegan", "ape", "mclust", "reshape2", "DESeq2", "ggplot2", "pheatmap",
  "RColorBrewer", "cluster", "dplyr", "PCAtools", "KOGMWU", "assertthat",
  "scales", "WGCNA", "tidyverse", "ggpubr", "rstatix", "MCMCglmm",
  "flashClust", "adegenet", "variancePartition", "BiocParallel", "limma", "car"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste(
      "Missing required R packages:",
      paste(missing_packages, collapse = ", "),
      "\nInstall them before running this script."
    )
  )
}

invisible(lapply(required_packages, library, character.only = TRUE))

# Global variables
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
message("Working directory set to: ", getwd())

# STEP 1 Data preparation & exploration ----

# coldata
coldata <- read.table(
  "data/metadata.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

rownames(coldata) <- coldata$Ind
coldata$Cross <- factor(coldata$Cross)

# import
cts <- read.table(
  "data/ortho_counts_matrix.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

rownames(cts) <- cts$ORTHO_ID
cts$ORTHO_ID <- NULL

dim(cts)        # rows = orthos (~20k), cols = samples (36)
head(cts[,1:4])

# force round
cts <- as.matrix(cts)
storage.mode(cts) <- "integer"
is.matrix(cts)            # TRUE
is.integer(cts)           # TRUE
any(cts < 0)              # FALSE
any(cts != floor(cts))    # FALSE"

any(cts < 0)                     # FALSE
any(cts != floor(cts))           # FALSE

cts<-cts[,colnames(cts) %in% rownames(coldata) ]
cts <- cts[,order(colnames(cts))]


## Add optionnal count
#build dds
dds <- DESeqDataSetFromMatrix(
  countData = cts,
  colData   = coldata,
  design    = ~ Cross
)


## filtering
dds <- estimateSizeFactors(dds)
dds <- estimateDispersions(dds)

idx <- rowSums(counts(dds,normalized=TRUE) >= 10 ) >= 3  # check use normalize
dds <- dds[idx,]
dim(dds)


## Estimate variables collinearity ##
cor(model.matrix(~ weight + Cross + gi, coldata))

X <- model.matrix(~ weight + Cross + gi, coldata)


# Convert matrix to dataframe (excluding intercept)
X_df <- as.data.frame(X[, -1])


# Fit an artificial regression model using the first variable as response
mod <- lm(as.formula(paste(colnames(X_df)[1], "~", paste(colnames(X_df)[-1], collapse = " + "))), data = X_df)

# Compute VIF
vif_values <- vif(mod)
print(vif_values)
#CrossAG  CrossGA  CrossGG       gi 
#1.485732 1.586618 2.223132 1.665512 
# explore metadata

ggplot(coldata,aes(x=Cross,y=gi,fill=Cross))+geom_boxplot() + geom_point(size=2,alpha=0.6)
ggplot(coldata,aes(x=Cross,y=weight,fill=Cross))+geom_boxplot() + geom_point(size=2,alpha=0.6)



## Check redundancy
majSequences <- function(counts, n=3, group, col=c("lightblue","orange","MediumVioletRed","SpringGreen"), outfile=TRUE){
  
  seqnames <- apply(counts, 2, function(x){x <- sort(x, decreasing=TRUE); names(x)[1:n]})
  seqnames <- unique(unlist(as.character(seqnames)))
  
  sum <- apply(counts,2,sum)
  counts <- counts[seqnames,]
  sum <- matrix(sum,nrow(counts),ncol(counts),byrow=TRUE)
  p <- round(100*counts/sum,digits=3)
  
  if (outfile) png(filename="output/check.maj.seq.ortho.png",width=min(3600,1800+800*ncol(counts)/10),height=1800,res=300)
  maj <- apply(p, 2, max)
  seqname <- rownames(p)[apply(p, 2, which.max)]
  x <- barplot(maj, col=col[as.integer(group)], main="Percentage of reads from most expressed sequence",
               ylim=c(0, max(maj)*1.2), las=2, ylab="Percentage of reads")
  legend("topright", levels(group), fill=col[1:nlevels(group)], bty="n")
  for (i in 1:length(seqname)) text(x[i], maj[i]/2, seqname[i], cex=0.8, srt=90, adj=0)
  if (outfile) dev.off()
  
  return(invisible(p))
}

majSequences(counts(dds), group=as.factor(coldata$Cross))

## Matrix log & VST transformations
### log norm
norm.counts <- counts(dds, normalized=TRUE)
log.norm.counts <- log2(norm.counts + 2)
df.log.norm.counts<-as.data.frame(log.norm.counts)
df.norm.counts<-as.data.frame(norm.counts)

### VST
vsd.fast <- vst(dds, fitType='local',blind=TRUE)
write.table(assay(vsd.fast),file="output/vst_matrix_ortho_2026.txt", quote=F)
vst.mat<-assay(vsd.fast)

## plot PCA & loadings
# plot options
theme_glob<-theme(axis.text.x=element_text(colour="black",size=12),
                  axis.text.y=element_text(colour="black",size=12),
                  axis.title.x=element_text(colour="black",size=12,face="bold"),
                  axis.title.y=element_text(colour="black",size=12,face="bold"),
                  panel.background = element_blank(),
                  panel.border=element_rect(colour = "black", fill=NA, size=1.5),
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(),
                  plot.margin=unit(c(1,1,1,1),"line"),
                  legend.title=element_blank(),
                  aspect.ratio=1,
                  legend.background = element_blank(),
                  legend.key = element_rect(fill=NA),
                  legend.text = element_text(colour="black",size=12,face="bold"))

couleurs=c("dodgerblue2","tomato2","orange2","black") 

size=c(2,3,5,6)

#plot global
plotData<-plotPCA(vsd.fast,intgroup=c("Cross","gi","weight"), returnData=TRUE)
percentVar <- round(100 * attr(plotData, "percentVar"))
PCA<-ggplot(plotData, aes(PC1, PC2,fill=Cross)) + geom_point(size=3, shape=21,alpha=0.7,stroke=1)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +coord_fixed()
PCA + geom_text(aes(label = name))

set.seed(1)
graph.pca<-PCA +theme_glob + scale_color_manual(values=couleurs)+
  scale_fill_manual(values = c("AA" = "tomato2", "GG" = "blue", "AG" = "orange2","GA"="grey50"),
                    labels = c("AA" = "Cangulata", "GG" = "Cgigas", "AG" = "Cang-Cgig","GA"="Cgig-Cang"))+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  annotate(geom="text", x=-10, y=40, label="PCA - orthogroups (20k entries 1:1)",
           color="black",size=4)
graph.pca
#dev.off()


# plot only nuclear genes only nuclear orthologs


mito <-read.table("data/mito.txt")
mito_ids <- mito$V1
keep_nuc <- !(rownames(vsd.fast) %in% mito_ids)

# PCA without mito
plotDataNuc<-plotPCA(vsd.fast[keep_nuc, ], intgroup=c("Cross","gi","weight"), returnData=TRUE)
percentVarNuc <- round(100 * attr(plotDataNuc, "percentVar"))
PCANuc<-ggplot(plotDataNuc, aes(PC1, PC2,fill=Cross)) + geom_point(size=3, shape=21,alpha=0.7,stroke=1)+
  xlab(paste0("PC1: ",percentVarNuc[1],"% variance")) +
  ylab(paste0("PC2: ",percentVarNuc[2],"% variance")) +coord_fixed()
PCANuc + geom_text(aes(label = name))

set.seed(1)
graph.pcaNuc<-PCANuc +theme_glob + scale_color_manual(values=couleurs)+
  scale_fill_manual(values = c("AA" = "tomato2", "GG" = "blue", "AG" = "orange2","GA"="grey50"),
                    labels = c("AA" = "Cangulata", "GG" = "Cgigas", "AG" = "Cang-Cgig","GA"="Cgig-Cang"))+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  annotate(geom="text", x=-10, y=30, label="PCA - Nucl. orthogroups (N = ~20k)",
           color="black",size=4)
graph.pcaNuc

# plot mito only
keep_mito <- rownames(vsd.fast) %in% mito_ids
plotDataMito<-plotPCA(vsd.fast[keep_mito, ], intgroup=c("Cross","gi","weight"), returnData=TRUE)
percentVarMito <- round(100 * attr(plotDataMito, "percentVar"))
PCAMito<-ggplot(plotDataMito, aes(PC1, PC2,fill=Cross)) + geom_point(size=3, shape=21,alpha=0.7,stroke=1)+
  xlab(paste0("PC1: ",percentVarMito[1],"% variance")) +
  ylab(paste0("PC2: ",percentVarMito[2],"% variance")) +coord_fixed()
PCAMito + geom_text(aes(label = name))

set.seed(1)
graph.pcaMito<-PCAMito +theme_glob + scale_color_manual(values=couleurs)+
  scale_fill_manual(values = c("AA" = "tomato2", "GG" = "blue", "AG" = "orange2","GA"="grey50"),
                    labels = c("AA" = "Cangulata", "GG" = "Cgigas", "AG" = "Cang-Cgig","GA"="Cgig-Cang"))+
  geom_hline(yintercept=0, linetype="dotted") +
  geom_vline(xintercept=0, linetype="dotted") +
  annotate(geom="text", x=-0, y=3, label="Mito. orthogroups (N = 37)",
           color="black",size=4)
graph.pcaMito


# Combine plots
library(patchwork)
combined_plot_fig1 <- (graph.pcaNuc | graph.pcaMito) +
  plot_annotation(tag_levels = "A")

combined_plot_fig1

ggsave(
  filename = "output/fig1_combined_pca.png",
  plot = combined_plot_fig1,
  width = 12,
  height = 6,
  dpi = 300
)

ggsave(
  filename = "output/fig1_combined_pca.pdf",
  plot = combined_plot_fig1,
  width = 12,
  height = 6
)

# plot Pheatmap of mitochrondria
# keep only mitochondrial genes
vsd_mito <- vst.mat[rownames(vst.mat) %in% mito_ids, ]
vsd_mito_scaled <- t(scale(t(vsd_mito)))

# extract sample metadata
coldata <- as.data.frame(colData(vsd.fast))

# example: genotype column called "cross"
# AG / AA -> mother = A
# GA / GG -> mother = G
coldata$mother <- ifelse(
  coldata$Cross %in% c("AA", "AG"),
  "A",
  "G"
)
coldata$mother <- factor(coldata$mother, levels = c("A", "G"))
# keep only annotation columns you want
table(coldata$Cross, coldata$mother)
annotation_col <- data.frame(
  mother = coldata$mother
)
rownames(annotation_col) <- rownames(coldata)

ann_colors <- list(
  mother = c(A = "#E69F00", G = "#56B4E9")
)

library(pheatmap)

pheatmap(
  vsd_mito_scaled,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  show_rownames = TRUE,
  fontsize_row = 6,
  fontsize_col = 10,
  main = "Mitochondrial genes (VST) – maternal origin"
)


## plot PCAtools
library(PCAtools)
# plot PCA
p <- PCAtools::pca(vst.mat, metadata = coldata, removeVar = 0.01)

biplot(p,showLoadings = TRUE,colby='Cross',
       vline=0, hline=0,
       legendPosition='right')

biplot(p,showLoadings = F,colby='Cross',
       vline=0, hline=0,
       legendPosition='right')

pairsplot(p,colby='Cross')
plotloadings(p)
screeplot(p)
ordered_loadings <- p$loadings[order(p$loadings[, "PC1"]), ]
head(ordered_loadings)

save(vst.mat,dds,coldata,mito_ids,file="output/step1.hybrids.ortho.Rda")



# STEP 2 WGCNA ----
ls()
rm(list=ls())
ls()

library(WGCNA)

load("output/step1.hybrids.ortho.Rda")

# Take a quick look at what is in the data set:

rm_samples <- grep("^(AA|GG)", colnames(dds), value = TRUE)
dds_hyb <- dds[!rownames(dds) %in% mito_ids, !colnames(dds) %in% rm_samples]
dds_hyb <- estimateSizeFactors(dds_hyb)
vsd.pure <- vst(dds_hyb, fitType='local',blind=TRUE)
vst.pure<-assay(vsd.pure)

datExpr0 = as.data.frame(t(vst.pure))
head((datExpr0))

# remove mito and keep only hybrid

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK



if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}


sampleTree = hclust(dist(datExpr0), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)


# Plot a line to show the cut
abline(h = 150, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 150, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
nSamples
variancedatExpr=as.vector(apply(as.matrix(datExpr),2,var, na.rm=T))
no.presentdatExpr=as.vector(apply(!is.na(as.matrix(datExpr)),2, sum) )
# Another way of summarizing the number of pressent entries
table(no.presentdatExpr)


# Keep only genes whose variance is non-zero and have at least 4 present entries
KeepGenes= variancedatExpr>0.05 # TEST for 0.05
table(KeepGenes)
datExpr=datExpr[, KeepGenes]

name_datExpr <-colnames(datExpr)

allTraits = read.table("data/metadata.txt",sep="\t", header=T,na.strings="NA");
names(allTraits)
summary(allTraits)
head(allTraits)
rownames(allTraits)<-allTraits$Ind
# Form a data frame analogous to expression data that will hold the clinical traits.
datExpr <- datExpr[ order(row.names(datExpr)), ]
allTraits <- allTraits[ order(row.names(allTraits)), ]
rownames(allTraits)

# refit column names
#Define desired order
new_order <- c("Ind", "Cross" , "gi","weight","AA","GG","AG","GA","matA","matG")

# Reorder columns
allTraits <- allTraits %>% select(all_of(new_order))
#allTraits$Ind<-NULL
allTraits$Cross<-as.factor(allTraits$Cross)
str(allTraits)
femaleSamples = rownames(datExpr);
traitRows = match(femaleSamples, allTraits$Ind);
datTraits = allTraits[traitRows, -1];
rownames(datTraits) = allTraits[traitRows, 1];
str(datTraits)
row.names(datTraits)

datTraits$Cross <- as.numeric(datTraits$Cross)
collectGarbage();

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits,signed= FALSE);
# Plot the sample dendrogram and the colors underneath.
#pdf("dendo_heatmap_subset.pdf",width=12,height=9)
par(mar=c(1, 10, 1, 1))
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

dev.off()

save(datExpr, datTraits, file = "output/dataInput_hybrids_ortho.Rda")

#setting is important, do not omit.
options(stringsAsFactors = FALSE);
# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
#enableWGCNAThreads()
# Load the data saved in the first part
lnames = load(file = "output/dataInput_hybrids_ortho.Rda");
#The variable lnames contains the names of loaded variables.
lnames

# build working matrix alias
datExpr_all <- datExpr

## Choose soft threshold power (recommended). If you already have a power, set it.
powers <- c(1:10, seq(12, 30, 2))
sft <- pickSoftThreshold(datExpr_all, powerVector = powers, verbose = 5)

# Plot strategy (same as WGCNA tutorial)
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 = 0.9
plot(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  xlab = "Soft Threshold (power)",
  ylab = "Scale Free Topology Model Fit, signed R^2",
  type = "n",
  main = "Scale independence"
)
text(
  sft$fitIndices[, 1],
  -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
  labels = powers,
  cex = cex1,
  col = "red"
)
abline(h = 0.9, col = "red")

plot(
  sft$fitIndices[, 1],
  sft$fitIndices[, 5],
  xlab = "Soft Threshold (power)",
  ylab = "Mean Connectivity",
  type = "n",
  main = "Mean connectivity"
)
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")

save(sft, file = "output/sft_hybrids_ortho.Rda")

# Pick a power; common: smallest with scale-free fit ~0.8
# If ambiguous, pick a reasonable value (e.g., 6–12 for signed networks).
#softPower <- if (!is.na(sft$powerEstimate)) sft$powerEstimate else 8

softPower = 10

cor <- WGCNA::cor
bicor <- WGCNA::bicor

net_all <- blockwiseModules(
  datExpr_all,
  power              = softPower,
  TOMType            = "signed",
  networkType        = "signed",
  minModuleSize      = 50,
  mergeCutHeight     = 0.25,
  numericLabels      = FALSE,
  pamRespectsDendro  = FALSE,
  saveTOMs           = FALSE,
  verbose            = 3,
  corType            = "pearson",
  corFnc             = WGCNA::cor,
  corOptions         = list(use = "p")
)

# net_all$colors are already colors because numericLabels = FALSE
moduleColors_all <- net_all$colors          # keep as-is
# no renaming
stopifnot(identical(moduleColors_all, net_all$colors))

length(net_all$colors)
ncol(datExpr_all)
is.null(names(net_all$colors))
head(names(net_all$colors))
head(colnames(datExpr_all))
moduleColors_all <- net_all$colors
names(moduleColors_all) <- colnames(datExpr_all)   # IMPORTANT: name by gene IDs
# Are they the same set of genes?
setequal(names(net_all$colors), colnames(datExpr_all))

# Do they match in the same order?
identical(names(net_all$colors), colnames(datExpr_all))

table(net_all$colors)[1:10]
table(moduleColors_all)[1:10]
stopifnot(all(net_all$colors == moduleColors_all))

MEs_all <- net_all$MEs
MEs_all <- orderMEs(MEs_all)

moduleColors <- moduleColors_all
MEs <- MEs_all
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder) - 1
geneTree <- net_all$dendrograms[[1]]
save(MEs, moduleLabels, moduleColors, geneTree, file = "output/networkConstruction_hybrids_sft10.Rda")


## 3) GLOBAL MODULE–TRAIT HEATMAP (Spearman, WGCNA-style)

## Convert phenotype and site to suitable trait matrix
coldata <- allTraits[rownames(datExpr_all), , drop = FALSE]

if (!"phenotype" %in% colnames(coldata)) {
  coldata$phenotype <- coldata$Cross
}
if (!"site" %in% colnames(coldata)) {
  coldata$site <- coldata$Cross
}
if (!"lustre" %in% colnames(coldata)) {
  coldata$lustre <- if ("gi" %in% colnames(coldata)) coldata$gi else NA_real_
}
if (!"quality" %in% colnames(coldata)) {
  coldata$quality <- if ("weight" %in% colnames(coldata)) coldata$weight else NA_real_
}
if (!"verif" %in% colnames(coldata)) {
  coldata$verif <- if ("matG" %in% colnames(coldata)) coldata$matG else 0
}

## Ensure factors
coldata$phenotype <- factor(coldata$phenotype)
coldata$site      <- factor(coldata$site)

## 1) phenotype -> single binary column (choose which level is 1)
lvl <- levels(coldata$phenotype)
if (length(lvl) < 2) {
  stop("`phenotype` needs at least 2 levels for module-trait modeling.")
}
# Set first level as reference; make second level = 1
coldata$phenotype01 <- as.numeric(coldata$phenotype == lvl[2])

## 2) site -> dummy columns (no intercept)
siteMat <- model.matrix(~ 0 + site, data = coldata)
colnames(siteMat) <- sub("^site", "site_", colnames(siteMat))

## 3) numeric / ordinal traits
coldata$lustre  <- suppressWarnings(as.numeric(coldata$lustre))
coldata$quality <- suppressWarnings(as.numeric(coldata$quality))
coldata$verif   <- suppressWarnings(as.numeric(coldata$verif))

## 4) combine into traitMat
traitMat <- cbind(
  phenotype = coldata$phenotype01,
  lustre    = coldata$lustre,
  quality   = coldata$quality,
  verif     = coldata$verif,
  siteMat
)

## 5) align with MEs
traitMat <- traitMat[rownames(MEs_all), , drop = FALSE]
stopifnot(identical(rownames(traitMat), rownames(MEs_all)))


moduleTraitCor <- cor(MEs_all, traitMat, method = "spearman", use = "pairwise.complete.obs")
moduleTraitP   <- corPvalueStudent(moduleTraitCor, nSamples = nrow(MEs_all))

textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(",
                     signif(moduleTraitP, 1), ")")
dim(textMatrix) <- dim(moduleTraitCor)

par(mar = c(6, 10, 3, 3))
labeledHeatmap(
  Matrix      = moduleTraitCor,
  xLabels     = colnames(traitMat),
  yLabels     = colnames(MEs_all),
  ySymbols    = colnames(MEs_all),
  colorLabels = FALSE,
  colors      = blueWhiteRed(50),
  textMatrix  = textMatrix,
  setStdMargins = FALSE,
  cex.text    = 0.7,
  zlim        = c(-1, 1),
  main        = "Global module–trait relationships (Spearman)"
)

## Optional: also compute an overall ANOVA p-value per module for phenotype (clean inference)
p_anova <- sapply(colnames(MEs_all), function(me) {
  fit <- lm(MEs_all[[me]] ~ coldata$phenotype)
  anova(fit)[1, "Pr(>F)"]
})
p_anova_fdr <- p.adjust(p_anova, method = "fdr")
anova_table <- data.frame(
  module = names(p_anova),
  p = p_anova,
  fdr = p_anova_fdr,
  stringsAsFactors = FALSE
)
anova_table <- anova_table[order(anova_table$fdr), ]
print(head(anova_table, 20))


## Plot with ggplot
library(tidyr)
library(ggplot2)
library(tibble)

# 1. Remove unwanted traits
traitMat_filt <- traitMat[, !colnames(traitMat) %in% c("lustre", "quality", "verif"), drop = FALSE]

# Keep only shared samples, and force identical order
common_samples <- intersect(rownames(MEs_all), rownames(traitMat_filt))

MEs_use   <- MEs_all[common_samples, , drop = FALSE]
trait_use <- traitMat_filt[common_samples, , drop = FALSE]

# Safety check
stopifnot(identical(rownames(MEs_use), rownames(trait_use)))
stopifnot(nrow(MEs_use) == nrow(trait_use))

# Correlation matrix
moduleTraitCor <- cor(
  MEs_use,
  trait_use,
  use = "pairwise.complete.obs",
  method = "spearman"
)

# P-value matrix
moduleTraitPvalue <- matrix(
  NA,
  nrow = ncol(MEs_use),
  ncol = ncol(trait_use),
  dimnames = list(colnames(MEs_use), colnames(trait_use))
)

for (i in seq_len(ncol(MEs_use))) {
  for (j in seq_len(ncol(trait_use))) {
    x <- MEs_use[, i]
    y <- trait_use[, j]

    ok <- complete.cases(x, y)

    if (sum(ok) >= 3) {
      tmp <- suppressWarnings(
        cor.test(x[ok], y[ok], method = "spearman", exact = FALSE)
      )
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


# 4. Convert to long format for ggplot
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

# Optional: order modules/traits
heatmap_df$module <- factor(heatmap_df$module, levels = rev(colnames(MEs_use)))
heatmap_df$trait  <- factor(heatmap_df$trait, levels = colnames(trait_use))
heatmap_df$module <- gsub("^ME", "", heatmap_df$module)

# 5. ggplot heatmap
p_heatmap <- ggplot(heatmap_df, aes(x = trait, y = module, fill = correlation)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 3, lineheight = 0.9) +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    limits = c(-1, 1),
    name = "Spearman\nrho"
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.ticks = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)
  )

p_heatmap

nGenes = ncol(datExpr)
nSamples = nrow(datExpr)
merge.ME <- merge(datTraits, MEs, by = 0)
write.table(merge.ME, file = "output/module_hybrids_cgigas.txt", quote = FALSE)



## Extract ME infos & annotations
#### Get info
cross = as.data.frame(datTraits$Cross);
names(cross) = "cross"
str(cross)
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
str(geneModuleMembership)


MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
str(MMPvalue)
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, cross, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(cross), sep="");
names(GSPvalue) = paste("p.GS.", names(cross), sep="");

## plot correlation
### significance by module

module = "red"
column = match(module, modNames);
moduleGenes = moduleColors==module;

pdf(file = "output/module_red.pdf", width = 6, height = 6)
plot<-verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, 1]),
                         xlab = paste("Module Membership in ", module, "module"),
                         ylab = "Gene significance for cross",
                         main = NULL, cex.lab = 1.2, cex.axis = 1.2, col = module)
plot + text(0.3, 0.78, "abs(cor) = 0.59\np-value < 0.001",cex = 1.5,font=2)
dev.off()



# export geneinfo
annot = read.table(file = "wgcna.annot.tab",sep="\t",header=T);
dim(annot)
names(annot)
probes = names(datExpr)
probes
probes2annot = match(probes, annot$Name)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

# Create the starting data frame
geneInfo0 = data.frame(geneID = probes,
                       geneSymbol = annot$Sp[probes2annot],
                       LocusLinkID = annot$ID[probes2annot],
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for temperature
modOrder = order(-abs(cor(MEs, cross, use = "p")));

# get relevant modules
# Get the corresponding Locuis Link IDs
allLLIDs = annot$name[probes2annot];
# $ Choose interesting modules
intModules = c("cyan","darkgrey","lightgreen","midnightblue","turquoise",
               "red","blue","brown","tan","darkgreen","black","darkred",
               "pink","grey60","yellow")
for (module in intModules)
{
  # Select module probes
  modGenes = (moduleColors==module)
  # Get their entrez ID codes
  modLLIDs = allLLIDs[modGenes];
  # Write them into a file
  fileName = paste("module_global_", module, ".txt", sep="");
  write.table(as.data.frame(modLLIDs), file = fileName,
              row.names = FALSE, col.names = FALSE,quote=FALSE)
}

# A
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.cross));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneinfo_hybrids_cgigas_cross.csv")

### Plot major
rownames(merge.ME)<-merge.ME$Row.names

mrg.good<-merge(coldata,merge.ME,by=0)

ggplot(mrg.good,aes(x=Cross.x,y=MEturquoise,fill=Cross.x))+
  geom_point(size=3,alpha=0.8,shape=21)

ggplot(mrg.good,aes(x=Cross.x,y=MEblue,fill=Cross.x))+
  geom_point(size=3,alpha=0.8,shape=21)

ggplot(mrg.good,aes(x=Cross.x,y=MEgreen,fill=gi.x))+
  geom_point(size=3,alpha=0.8,shape=21)

ggplot(merge.ME,aes(x=as.factor(condition),y=MEdarkorange2,color=as.factor(condition),shape=as.factor(tank)))+
  geom_point(size=4)

## STEP 3 DGE a,d dominance effect ----

load("step1.hybrids.ortho.Rda")

# Nuclear
dds_hyb <- dds[!rownames(dds) %in% mito_ids, ]
### test standard
design(dds_hyb) <- ~ Cross
dds_hyb <- DESeq(dds_hyb)
resultsNames(dds_hyb)

# diff 
res_GG_AA <- results(dds_hyb, contrast = c("Cross","GG","AA"))
summary(res_GG_AA$log2FoldChange)
mean(abs(res_GG_AA$log2FoldChange) > 1, na.rm=TRUE) 
res_AG_AA <- results(dds_hyb, contrast = c("Cross","AG","AA"))
res_GA_AA <- results(dds_hyb, contrast = c("Cross","GA","AA"))


# compute dominance
D_AG <- (res_AG_AA$log2FoldChange - 0.5 * res_GG_AA$log2FoldChange) /
  (0.5 * res_GG_AA$log2FoldChange)

D_GA <- (res_GA_AA$log2FoldChange - 0.5 * res_GG_AA$log2FoldChange) /
  (0.5 * res_GG_AA$log2FoldChange)

# interpreation

#  D	Pattern
#~0	Additive
#~+1	GG-dominant
#~−1	AA-dominant
#> +1	Overdominant
#< −1	Underdominant

# Explore in AG
alpha <- 0.05

class_AG <- rep("Ambiguous", nrow(res_AG_AA))

# Additive
class_AG[
  res_AG_AA$padj < alpha &
    res_GG_AA$padj < alpha &
    abs(res_AG_AA$log2FoldChange - 0.5 * res_GG_AA$log2FoldChange) < 0.25
] <- "Additive"

# GG-dominant
class_AG[
  res_AG_AA$padj < alpha &
    abs(res_AG_AA$log2FoldChange - res_GG_AA$log2FoldChange) < 0.25
] <- "GG_dominant"

# AA-dominant
class_AG[
  res_AG_AA$padj < alpha &
    abs(res_AG_AA$log2FoldChange) < 0.25
] <- "AA_dominant"

# Overdominant
class_AG[
  res_AG_AA$padj < alpha &
    res_GG_AA$padj < alpha &
    res_AG_AA$log2FoldChange > res_GG_AA$log2FoldChange + 0.25
] <- "Overdominant"

# Underdominant
class_AG[
  res_AG_AA$padj < alpha &
    res_GG_AA$padj < alpha &
    res_AG_AA$log2FoldChange < -0.25
] <- "Underdominant"


table(class_AG)
prop.table(table(class_AG))
# same for GA
alpha <- 0.05

class_GA <- rep("Ambiguous", nrow(res_GA_AA))

# Additive (mid-parent)
class_GA[
  res_GA_AA$padj < alpha &
    res_GG_AA$padj < alpha &
    abs(res_GA_AA$log2FoldChange - 0.5 * res_GG_AA$log2FoldChange) < 0.25
] <- "Additive"

# GG-dominant
class_GA[
  res_GA_AA$padj < alpha &
    abs(res_GA_AA$log2FoldChange - res_GG_AA$log2FoldChange) < 0.25
] <- "GG_dominant"

# AA-dominant
class_GA[
  res_GA_AA$padj < alpha &
    abs(res_GA_AA$log2FoldChange) < 0.25
] <- "AA_dominant"

# Overdominant
class_GA[
  res_GA_AA$padj < alpha &
    res_GG_AA$padj < alpha &
    res_GA_AA$log2FoldChange > res_GG_AA$log2FoldChange + 0.25
] <- "Overdominant"

# Underdominant
class_GA[
  res_GA_AA$padj < alpha &
    res_GG_AA$padj < alpha &
    res_GA_AA$log2FoldChange < -0.25
] <- "Underdominant"

# Summary
table(class_GA)


# compare tables
table(class_AG, class_GA)
# results
#Ambiguous → GG_dominant     = 50
#Ambiguous → Overdominant    = 15
#Ambiguous → Underdominant   = 95
# interpretations:
#Maternal effects/ Mito–nuclear interactions / Trans-regulatory asymmetr
#Dominance patterns are largely conserved between reciprocal hybrids, 
#yet ~1–2% of expressed genes exhibit directional dominance or misexpression, suggesting parent-of-origin–dependent regulatory divergence.”

# explore asymetricale genes
asymmetric <- class_AG != class_GA & 
  !(class_AG == "Ambiguous" & class_GA == "Ambiguous")

sum(asymmetric)
mean(asymmetric)
table(class_AG[asymmetric])
table(class_GA[asymmetric])

# plot pheatmap
library(pheatmap)

mat <- as.matrix(table(class_AG, class_GA))
pheatmap(mat / rowSums(mat),
         display_numbers = mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE)


## plot dominance
plot(res_GG_AA$log2FoldChange,
     res_AG_AA$log2FoldChange,
     pch=16, cex=0.4,
     xlab="GG vs AA",
     ylab="AG vs AA")
abline(0,1,col="red")     # GG dominance
abline(0,0.5,col="blue") # additive
abline(0,0,col="green")  # AA dominance

# classify genes to the distance of the line
delta <- 0.25

class <- rep("Other", nrow(res_AG_AA))

class[abs(res_AG_AA$log2FoldChange - 0.5 * res_GG_AA$log2FoldChange) < delta] <- "Additive"
class[abs(res_AG_AA$log2FoldChange - res_GG_AA$log2FoldChange) < delta] <- "GG_dominant"
class[abs(res_AG_AA$log2FoldChange) < delta] <- "AA_dominant"

class[res_AG_AA$log2FoldChange > res_GG_AA$log2FoldChange + delta] <- "Overdominant"
class[res_AG_AA$log2FoldChange < -delta] <- "Underdominant"

## plot individul gene
data.gene.merged<-as.data.frame(merge(t(vst.mat),coldata,by=0))

#sortilin-like
ggplot(data.gene.merged,aes(x=Cross,y=ORTHO_015690,fill=as.factor(Cross)))+
  geom_point(size=3,shape=21,alpha=0.6) 

#translocation protein SEC63
ggplot(data.gene.merged,aes(x=Cross,y=`gene-LOC128169801`,fill=as.factor(Cross)))+
  geom_point(size=3,shape=21,alpha=0.6)


## For mito only
# Nuclear
dds_mit <- dds[rownames(dds) %in% mito_ids, ]
dim(dds_mit)
### test standard
design(dds_mit) <- ~ Cross
dds_mit <- DESeq(dds_mit)
resultsNames(dds_mit)

# diff 
resmit_GG_AA <- results(dds_mit, contrast = c("Cross","GG","AA"))
resmit_AG_AA <- results(dds_mit, contrast = c("Cross","AG","AA"))
resmit_GA_AA <- results(dds_mit, contrast = c("Cross","GA","AA"))


# compute dominance
D_AG <- (resmit_AG_AA$log2FoldChange - 0.5 * resmit_GG_AA$log2FoldChange) /
  (0.5 * resmit_GG_AA$log2FoldChange)

D_GA <- (resmit_GA_AA$log2FoldChange - 0.5 * resmit_GG_AA$log2FoldChange) /
  (0.5 * resmit_GG_AA$log2FoldChange)

# interpreation

#  D	Pattern
#~0	Additive
#~+1	GG-dominant
#~−1	AA-dominant
#> +1	Overdominant
#< −1	Underdominant

# Explore in AG
alpha <- 0.05

class_AG <- rep("Ambiguous", nrow(resmit_AG_AA))

# Additive
class_AG[
  resmit_AG_AA$padj < alpha &
    resmit_GG_AA$padj < alpha &
    abs(resmit_AG_AA$log2FoldChange - 0.5 * resmit_GG_AA$log2FoldChange) < 0.25
] <- "Additive"

# GG-dominant
class_AG[
  resmit_AG_AA$padj < alpha &
    abs(resmit_AG_AA$log2FoldChange - resmit_GG_AA$log2FoldChange) < 0.25
] <- "GG_dominant"

# AA-dominant
class_AG[
  resmit_AG_AA$padj < alpha &
    abs(resmit_AG_AA$log2FoldChange) < 0.25
] <- "AA_dominant"

# Overdominant
class_AG[
  resmit_AG_AA$padj < alpha &
    resmit_GG_AA$padj < alpha &
    resmit_AG_AA$log2FoldChange > resmit_GG_AA$log2FoldChange + 0.25
] <- "Overdominant"

# Underdominant
class_AG[
  resmit_AG_AA$padj < alpha &
    resmit_GG_AA$padj < alpha &
    resmit_AG_AA$log2FoldChange < -0.25
] <- "Underdominant"


table(class_AG)
prop.table(table(class_AG))
# same for GA
alpha <- 0.05

class_GA <- rep("Ambiguous", nrow(resmit_GA_AA))

# Additive (mid-parent)
class_GA[
  resmit_GA_AA$padj < alpha &
    resmit_GG_AA$padj < alpha &
    abs(resmit_GA_AA$log2FoldChange - 0.5 * resmit_GG_AA$log2FoldChange) < 0.25
] <- "Additive"

# GG-dominant
class_GA[
  resmit_GA_AA$padj < alpha &
    abs(resmit_GA_AA$log2FoldChange - resmit_GG_AA$log2FoldChange) < 0.25
] <- "GG_dominant"

# AA-dominant
class_GA[
  resmit_GA_AA$padj < alpha &
    abs(resmit_GA_AA$log2FoldChange) < 0.25
] <- "AA_dominant"

# Overdominant
class_GA[
  resmit_GA_AA$padj < alpha &
    resmit_GG_AA$padj < alpha &
    resmit_GA_AA$log2FoldChange > resmit_GG_AA$log2FoldChange + 0.25
] <- "Overdominant"

# Underdominant
class_GA[
  resmit_GA_AA$padj < alpha &
    resmit_GG_AA$padj < alpha &
    resmit_GA_AA$log2FoldChange < -0.25
] <- "Underdominant"

# Summary
table(class_GA)


# compare tables
table(class_AG, class_GA)
# results
#Ambiguous → GG_dominant     = 50
#Ambiguous → Overdominant    = 15
#Ambiguous → Underdominant   = 95
# interpretations:
#Maternal effects/ Mito–nuclear interactions / Trans-regulatory asymmetr
#Dominance patterns are largely conserved between reciprocal hybrids, 
#yet ~1–2% of expressed genes exhibit directional dominance or misexpression, suggesting parent-of-origin–dependent regulatory divergence.”

# explore asymetricale genes
asymmetric <- class_AG != class_GA & 
  !(class_AG == "Ambiguous" & class_GA == "Ambiguous")

sum(asymmetric)
mean(asymmetric)
table(class_AG[asymmetric])
table(class_GA[asymmetric])

# plot pheatmap
library(pheatmap)

mat <- as.matrix(table(class_AG, class_GA))
pheatmap(mat / rowSums(mat),
         display_numbers = mat,
         cluster_rows = FALSE,
         cluster_cols = FALSE)


## plot dominance
plot(res_GG_AA$log2FoldChange,
     res_AG_AA$log2FoldChange,
     pch=16, cex=0.4,
     xlab="GG vs AA",
     ylab="AG vs AA")
abline(0,1,col="red")     # GG dominance
abline(0,0.5,col="blue") # additive
abline(0,0,col="green")  # AA dominance

# classify genes to the distance of the line
delta <- 0.25

class <- rep("Other", nrow(res_AG_AA))

class[abs(res_AG_AA$log2FoldChange - 0.5 * res_GG_AA$log2FoldChange) < delta] <- "Additive"
class[abs(res_AG_AA$log2FoldChange - res_GG_AA$log2FoldChange) < delta] <- "GG_dominant"
class[abs(res_AG_AA$log2FoldChange) < delta] <- "AA_dominant"

class[res_AG_AA$log2FoldChange > res_GG_AA$log2FoldChange + delta] <- "Overdominant"
class[res_AG_AA$log2FoldChange < -delta] <- "Underdominant"

## plot individul gene
data.gene.merged<-as.data.frame(merge(t(vst.mat),coldata,by=0))

#sortilin-like
ggplot(data.gene.merged,aes(x=Cross,y=ORTHO_015690,fill=as.factor(Cross)))+
  geom_point(size=3,shape=21,alpha=0.6) 

#translocation protein SEC63
ggplot(data.gene.merged,aes(x=Cross,y=`gene-LOC128169801`,fill=as.factor(Cross)))+
  geom_point(size=3,shape=21,alpha=0.6)

## cluster DEGS ----


# Log-transform RNA-seq counts for better Gaussian distribution

# get only deg in LRT
vst.mat.deg<-as.data.frame(vst.mat[rownames(vst.mat) %in% rownames(res.filt),])
#vst.mat.deg$Gene <- rownames(vst.mat.deg)  # Add gene names as a column

# Fit model-based clustering
mclust_result <- Mclust(vst.mat.deg)

# View clustering summary
summary(mclust_result)

# Extract cluster assignments
clusters <- mclust_result$classification

# Standardize gene expression (z-score transformation)
standardized_data <- t(scale(t(vst.mat.deg)))  # Scale genes (rows)

# Convert standardized data into a data frame
df <- as.data.frame(standardized_data)
df$Cluster <- as.factor(clusters)  # Add cluster information
df$GeneID <- rownames(df)  # Assign gene names

sample_map <- data.frame(Sample = rownames(coldata), Cross = coldata$Cross, row.names = NULL)

# Convert standardized gene expression into long format
df_long <- melt(df, id.vars = c("GeneID", "Cluster"), variable.name = "Sample", value.name = "Expression")

# Merge the treatment mapping based on the sample names
df_long <- left_join(df_long, sample_map, by = "Sample")

#summarize by treatment
df_cross <- df_long %>%
  group_by(GeneID, Cluster, Cross) %>%
  summarise(Expression = mean(Expression), .groups = "drop")

#compute eigengen value
eigengene_df <- df_cross %>%
  group_by(Cluster, Cross) %>%
  summarise(Eigengene = mean(Expression), .groups = "drop")


#plot
ggplot(df_cross, aes(x = Cross, y = Expression, group = GeneID)) +
  geom_line(color = "grey", alpha = 0.8) +  # Grey lines for all genes
  geom_line(data = eigengene_df, aes(x = Cross, y = Eigengene, group = Cluster, color = Cluster), size = 1) +  # Bold cluster eigengene
  facet_wrap(~Cluster, scales = "free_x") +  # Separate clusters into panels
  theme_minimal() +
  labs(title = "Parallel Coordinate Plot of Standardized Gene Expression",
       x = "Treatment",
       y = "Standardized Expression (Z-score)") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate x-axis labels

# Run PCA on VST data
pca_result <- PCA(t(standardized_data), scale.unit = TRUE, ncp = 10)  # Reduce to 10 PCs

# Extract first 2 PCs
pca_data <- as.data.frame(pca_result$ind$coord)
colnames(pca_data)[1:2] <- c("PC1", "PC2")

# Plot PCA to check cluster structure
ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point() +
  theme_minimal() +
  ggtitle("PCA of Gene Expression Data")



# STEP 4 : compute ePST values ----


vsd.fast <- vst(dds_hyb, fitType='local',blind=TRUE)
meta <- as.data.frame(colData(vsd.fast))

expr <- assay(vsd.fast)   # matrix: genes x samples
dim(expr)
keep <- meta$Cross %in% c("AA","GG")
expr2 <- expr[, keep, drop = FALSE]   # IMPORTANT drop=FALSE
meta2 <- meta[keep, , drop = FALSE]

dim(expr2)
table(meta2$Cross)

# define group
group <- factor(meta2$Cross)

# function ePST
compute_ePST <- function(x, group) {
  # x = expression values for one gene
  # group = factor with two levels (AA, GG)
  
  means <- tapply(x, group, mean)
  vars  <- tapply(x, group, var)
  
  V_between <- var(means)
  V_within  <- mean(vars, na.rm = TRUE)
  
  ePST <- V_between / (V_between + 2 * V_within)
  return(ePST)
}


ePST <- apply(expr2, 1, compute_ePST, group = group)
summary(ePST)
mean(ePST, na.rm = TRUE)
sd(ePST, na.rm = TRUE)

# ePST range | Interpretation                       |
#   | ---------- | ------------------------------------ |
#   | < 0.05     | Mostly plastic / shared regulation   |
#   | 0.05–0.12  | Weak differentiation                 |
#   | > 0.12     | Regulatory divergence                |
#   | > 0.3      | Strong, possibly selected divergence |

# outliers ePST
thr95 <- quantile(ePST, 0.95, na.rm=TRUE)
thr99 <- quantile(ePST, 0.99, na.rm=TRUE)

thr95; thr99
sum(ePST > thr95, na.rm=TRUE)

table(group)
mean(ePST > 0.12, na.rm=TRUE)
quantile(ePST, c(0.5, 0.9, 0.95, 0.99), na.rm=TRUE)

boxplot(ePST ~ class_AG)

# select genes of interest
## Strong divergent
set_strong <- res_GG_AA$padj < 1e-4 &
  abs(res_GG_AA$log2FoldChange) > 2 &
  ePST > quantile(ePST, 0.95)
set_strong
## consistent, but more subtile
set_consistent <- ePST > quantile(ePST, 0.95) &
  abs(res_GG_AA$log2FoldChange) < 1
## large shift but noisy
set_noisy <- abs(res_GG_AA$log2FoldChange) > 2 &
  ePST < quantile(ePST, 0.5)

# exp diff vs ePST
plot(abs(res_GG_AA$log2FoldChange), ePST,
     pch=16, cex=0.4,
     xlab="|log2FC| (GG vs AA)",
     ylab="ePST")
cor.test(abs(res_GG_AA$log2FoldChange), ePST, use="complete.obs", method="spearman")

## Overlay ePST and dominance
res <- as.data.frame(res_GG_AA)
res$ORTHO_ID <- rownames(res)

# add ePST
res$ePST <- ePST[res$ORTHO_ID]

# add class_AG (robust to either named or positional vector)
if (!is.null(names(class_AG))) {
  res$class_AG <- class_AG[res$ORTHO_ID]
} else {
  res$class_AG <- class_AG
}

# keep only complete cases for plot
plotdf <- res[complete.cases(res$log2FoldChange, res$ePST, res$class_AG), ]
plotdf$absLFC <- abs(plotdf$log2FoldChange)
plotdf$class_AG <- factor(plotdf$class_AG)

plotdf$class_AG <- factor(plotdf$class_AG,
                          levels=c("Ambiguous","Additive","AA_dominant","GG_dominant",
                                   "Overdominant","Underdominant"))

library(ggplot2)

ggplot(plotdf, aes(x = absLFC, y = ePST, color = class_AG)) +
  geom_point(size = 0.6, alpha = 0.8) +
  labs(x = "|log2FC| (GG vs AA)", y = "ePST", color = "AG class") +
  theme_classic()
thr95 <- quantile(plotdf$ePST, 0.95, na.rm=TRUE)

ggplot(plotdf, aes(absLFC, ePST, color=class_AG)) +
  geom_point(size=0.6, alpha=0.35) +
  geom_hline(yintercept = thr95, linetype="dashed") +
  geom_vline(xintercept = 2, linetype="dashed") +
  labs(x="|log2FC| (GG vs AA)", y="ePST") +
  theme_classic()

plotdf$hit <- with(plotdf, ePST > quantile(ePST, 0.95, na.rm=TRUE) & absLFC > 2 & padj < 1e-4)

ggplot(plotdf, aes(absLFC, ePST)) +
  geom_point(aes(color=class_AG), size=0.5, alpha=0.25) +
  geom_point(data=subset(plotdf, hit), size=1.2) +
  geom_hline(yintercept = thr95, linetype="dashed") +
  geom_vline(xintercept = 2, linetype="dashed") +
  labs(x="|log2FC| (GG vs AA)", y="ePST", color="AG class") +
  theme_classic()

ggplot(plotdf, aes(absLFC, ePST)) +
  geom_bin2d(bins = 80) +
  facet_wrap(~ class_AG) +
  labs(x="|log2FC| (GG vs AA)", y="ePST") +
  theme_classic()

# STEP 5 GOMWU ----

# Edit these to match your data file names: 
input="deg.fullmod.csv" # two columns of comma-separated values: gene id, continuous measure of significance. To perform standard GO enrichment analysis based on Fisher's exact test, use binary measure (0 or 1, i.e., either sgnificant or not).
goAnnotations="dlabrax.annotations.go.tab" # two-column, tab-delimited, one line per gene, multiple GO terms separated by semicolon. If you have multiple lines per gene, use nrify_GOtable.pl prior to running this script.
goDatabase="go.obo" # download 24Nov2023 from http://www.geneontology.org/GO.downloads.ontology.shtml
source("gomwu.functions.R")

# GO Biological processes
goDivision="MF" # either MF, or BP, or CC

# Calculating stats. It might take ~3 min for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="~/miniconda3/envs/perl/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=3,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #Alternative="g", # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           #Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.


# Plotting results


results=gomwuPlot(input,goAnnotations,goDivision,
                  #	absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  absValue=0.001,
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=1, # height of the hierarchical clustering tree
                  colors=c("lightgrey","grey","black") # these are default colors, un-remar and change if needed
)

# text representation of results, with actual adjusted p-values
sig.mwus=results[[1]]
resFileName <- paste(goDivision,"_",input,".txt",sep="")
write.table(sig.mwus,file=resFileName,quote=F)

# ------- extracting representative GOs

# this module chooses GO terms that best represent *independent* groups of significant GO terms

pcut=1e-2 # adjusted pvalue cutoff for representative GO
hcut=0.9 # height at which cut the GO terms tree to get "independent groups". 

# plotting the GO tree with the cut level (un-remark the next two lines to plot)
plot(results[[2]],cex=0.6)
abline(h=hcut,col="red")

# cutting
ct=cutree(results[[2]],h=hcut)
annots=c();ci=1
for (ci in unique(ct)) {
  message(ci)
  rn=names(ct)[ct==ci]
  obs=grep("obsolete",rn)
  if(length(obs)>0) { rn=rn[-obs] }
  if (length(rn)==0) {next}
  rr=results[[1]][rn,]
  bestrr=rr[which(rr$pval==min(rr$pval)),]
  best=1
  if(nrow(bestrr)>1) {
    nns=sub(" .+","",row.names(bestrr))
    fr=c()
    for (i in 1:length(nns)) { fr=c(fr,eval(parse(text=nns[i]))) }
    best=which(fr==max(fr))
  }
  if (bestrr$pval[best]<=pcut) { annots=c(annots,sub("\\d+\\/\\d+ ","",row.names(bestrr)[best]))}
}

mwus=read.table(paste("MWU",goDivision,input,sep="_"),header=T)

bestGOs=mwus[mwus$name %in% annots,]
bestGOs

csvFileName <- paste(goDivision,"_",input,".Rda",sep="")
save(results,file=csvFileName)


# GO cellular component
goDivision="BP" # either MF, or BP, or CC


# Calculating stats. It might take ~3 min for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="~/miniconda3/envs/perl/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=3,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #Alternative="g", # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.


# Plotting results


results=gomwuPlot(input,goAnnotations,goDivision,
                  #	absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  absValue=0.001,
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=1, # height of the hierarchical clustering tree
                  colors=c("lightgrey","grey","black") # these are default colors, un-remar and change if needed
)

# text representation of results, with actual adjusted p-values
sig.mwus=results[[1]]
resFileName <- paste(goDivision,"_",input,".txt",sep="")
write.table(sig.mwus,file=resFileName,quote=F)

# GO cellular component
goDivision="CC" # either MF, or BP, or CC


# Calculating stats. It might take ~3 min for MF and BP. Do not rerun it if you just want to replot the data with different cutoffs, go straight to gomwuPlot. If you change any of the numeric values below, delete the files that were generated in previos runs first.
gomwuStats(input, goDatabase, goAnnotations, goDivision,
           perlPath="~/miniconda3/envs/perl/bin/perl", # replace with full path to perl executable if it is not in your system's PATH already
           largest=0.1,  # a GO category will not be considered if it contains more than this fraction of the total number of genes
           smallest=3,   # a GO category should contain at least this many genes to be considered
           clusterCutHeight=0.25, # threshold for merging similar (gene-sharing) terms. See README for details.
           #Alternative="g", # by default the MWU test is two-tailed; specify "g" or "l" of you want to test for "greater" or "less" instead. 
           Module=TRUE,Alternative="g" # un-remark this if you are analyzing a SIGNED WGCNA module (values: 0 for not in module genes, kME for in-module genes). In the call to gomwuPlot below, specify absValue=0.001 (count number of "good genes" that fall into the module)
           #Module=TRUE # un-remark this if you are analyzing an UNSIGNED WGCNA module 
)
# do not continue if the printout shows that no GO terms pass 10% FDR.


# Plotting results


results=gomwuPlot(input,goAnnotations,goDivision,
                  #	absValue=-log(0.05,10),  # genes with the measure value exceeding this will be counted as "good genes". Specify absValue=0.001 if you are doing Fisher's exact test for standard GO enrichment or analyzing a WGCNA module (all non-zero genes = "good genes").
                  absValue=0.001,
                  level1=0.1, # FDR threshold for plotting. Specify level1=1 to plot all GO categories containing genes exceeding the absValue.
                  level2=0.05, # FDR cutoff to print in regular (not italic) font.
                  level3=0.01, # FDR cutoff to print in large bold font.
                  txtsize=1.2,    # decrease to fit more on one page, or increase (after rescaling the plot so the tree fits the text) for better "word cloud" effect
                  treeHeight=1, # height of the hierarchical clustering tree
                  colors=c("lightgrey","grey","black") # these are default colors, un-remar and change if needed
)

# text representation of results, with actual adjusted p-values
sig.mwus=results[[1]]
resFileName <- paste(goDivision,"_",input,".txt",sep="")
write.table(sig.mwus,file=resFileName,quote=F)

# this module chooses GO terms that best represent *independent* groups of significant GO terms

pcut=1e-2 # adjusted pvalue cutoff for representative GO
hcut=0.9 # height at which cut the GO terms tree to get "independent groups". 

# plotting the GO tree with the cut level (un-remark the next two lines to plot)
# plot(results[[2]],cex=0.6)
# abline(h=hcut,col="red")

# cutting
ct=cutree(results[[2]],h=hcut)
annots=c();ci=1
for (ci in unique(ct)) {
  message(ci)
  rn=names(ct)[ct==ci]
  obs=grep("obsolete",rn)
  if(length(obs)>0) { rn=rn[-obs] }
  if (length(rn)==0) {next}
  rr=results[[1]][rn,]
  bestrr=rr[which(rr$pval==min(rr$pval)),]
  best=1
  if(nrow(bestrr)>1) {
    nns=sub(" .+","",row.names(bestrr))
    fr=c()
    for (i in 1:length(nns)) { fr=c(fr,eval(parse(text=nns[i]))) }
    best=which(fr==max(fr))
  }
  if (bestrr$pval[best]<=pcut) { annots=c(annots,sub("\\d+\\/\\d+ ","",row.names(bestrr)[best]))}
}

mwus=read.table(paste("MWU",goDivision,input,sep="_"),header=T)
bestGOs=mwus[mwus$name %in% annots,]
bestGOs

# STEP 6 : Run RDA ----
library(vegan)
library(PCPS)
library(dplyr)

load("step1.hybrids.ortho.Rda")
genet=assay(vsd.fast)
genet_trans=t(genet)
coldata<-coldata[,c("Cross","gi","weight")]
merge.data.ordered <- merge(coldata,genet_trans,by=0)
merge.data.ordered$Row.names <- NULL

merge.data.ordered <- na.omit(merge.data.ordered)
genet.pcoa=pcoa(daisy(merge.data.ordered[,4:ncol(merge.data.ordered)], metric="euclidean"))
genet.pcoa$values
pcoa.sig(genet_trans, method = "gower")

subset.variables <- merge.data.ordered %>% select(Cross,gi,weight)

subset.variables$Cross<-as.numeric(as.factor(subset.variables$Cross))
subset.variables$gi <- as.numeric(subset.variables$gi)
subset.variables$weight <- as.numeric(subset.variables$weight)

subset.variables <- na.omit(subset.variables)

sapply(subset.variables, class)

# run modl
### Ordistep function allows to add and remove variable to maximise the explained variance. 
#To avoid overfitting, selected variables should not explained more than the global model (36%) 
#OrdiR2step will start working from an empty model without explanatory variables, just the intercept
rda0<-vegan::rda(genet.pcoa$vectors[,1:11] ~ 1, subset.variables)
#OrdiR2step will move towards the global model with all explanatory variables
rdaG<- vegan::rda(genet.pcoa$vectors[,1:11] ~ ., subset.variables)

set.seed(10); Sel <- ordiR2step(rda0, scope = formula(rdaG), direction="both")
Sel$anova
rda_genet=rda(formula = genet.pcoa$vectors[, 1:11] ~ ., subset.variables, scale = F)
anova(rda_genet, step=1000) #anova.cca(rda_genet)
anova(rda_genet, step=1000, by='margin') #by='terms'
anova(rda_genet, step=1000, by='axis') #Root(vif) >2 coll. is high

# check collineearity
sqrt(vif.cca(rda_genet))

RsquareAdj(rda_genet)

# check partial
rda_partial=rda(genet.pcoa$vectors[,1:11],subset.variables$Cross,cbind(subset.variables$gi + subset.variables$weight), scale=F)
rda_partial
anova(rda_partial, step=1000)
RsquareAdj(rda_partial)

rda_partial2=rda(genet.pcoa$vectors[,1:10],subset.variables$gi,cbind(subset.variables$Cross+ subset.variables$weight), scale=F)
rda_partial2
anova(rda_partial2, step=1000)
RsquareAdj(rda_partial2)

summary(coldata[,c(1,2)])


# STEP 7: Explore reads-origin (ASE) ----
library(data.table)

sp1 <- fread("freq_sp1.frq")
sp2 <- fread("freq_sp2.frq")

# Merge by SNP ID or position
merged <- merge(sp1, sp2, by = c("CHROM", "POS"))

# Parse allele frequencies
parse_freq <- function(freq_col) {
  as.numeric(sub(".*:", "", freq_col))
}

merged$AF1 <- parse_freq(merged$MAF.x)
merged$AF2 <- parse_freq(merged$MAF.y)

# Get SNPs with fixed or nearly fixed differences
diagnostic_snps <- merged[abs(AF1 - AF2) >= 0.9]


# STEP 8 Module preservation ----
# exploration based in ortholog matrix
ls()
rm(list=ls())
ls()

suppressPackageStartupMessages({
  library(WGCNA)
  library(dynamicTreeCut)
  library(netrep)
  library(dplyr)
})

setwd(project_dir)

load("output/step1.hybrids.ortho.Rda")

# Take a quick look at what is in the data set:

# Test WGCNA sft for AG 
rm_samples <- grep("^(AA|GG)", colnames(dds), value = TRUE)
dds_AG <- dds[!rownames(dds) %in% mito_ids, !colnames(dds) %in% rm_samples]
dds_AG <- estimateSizeFactors(dds_AG)
vsd.hyb <- vst(dds_AG, fitType='local',blind=TRUE)
vst.hyb<-assay(vsd.hyb)





options(stringsAsFactors = FALSE)
allowWGCNAThreads()  # or disable if on cluster with strict threading rules

# -----------------------------
# INPUTS YOU MUST PROVIDE
# -----------------------------
# vst: matrix of VST values, rows = genes, cols = samples
# meta: data.frame with rownames(meta) = sample IDs matching colnames(vst)
#       and a column 'cross' with values "AG" or "GA"

# Example checks:
stopifnot(all(colnames(vst.hyb) %in% rownames(coldata)))
coldata <- coldata[colnames(vst.hyb), , drop = FALSE]
stopifnot(all(coldata$Cross %in% c("AG","GA")))

# -----------------------------
# 0) Basic QC: remove bad samples/genes if needed
# -----------------------------
gsg <- goodSamplesGenes(t(vst.hyb), verbose = 3)  # WGCNA expects samples x genes, so transpose
if (!gsg$allOK) {
  message("Removing bad samples/genes flagged by goodSamplesGenes()")
  if (sum(!gsg$goodGenes) > 0) vst.hyb <- vst.hyb[gsg$goodGenes, ]
  if (sum(!gsg$goodSamples) > 0) {
    badSamps <- colnames(vst.hyb)[!gsg$goodSamples]
    vst.hyb <- vst.hyb[, gsg$goodSamples]
    coldata <- coldata[colnames(vst.hyb), , drop = FALSE]
    message("Dropped samples: ", paste(badSamps, collapse = ", "))
  }
}

# -----------------------------
# 1) Choose a shared gene set (critical for netrep)
#    With n=9 per group, DO NOT use all genes.
# -----------------------------
# Option A: top variable genes across ALL samples
topN <- 8000  # adjust (3000–8000 is a reasonable range for n=9/group)
gene_var <- apply(vst.hyb, 1, var, na.rm = TRUE)
keep_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(min(topN, length(gene_var)))]
vst_f <- vst.hyb[keep_genes, ]

# -----------------------------
# 2) Split by cross
# -----------------------------
vst_AG <- vst_f[, coldata$Cross == "AG", drop = FALSE]
vst_GA <- vst_f[, coldata$Cross == "GA", drop = FALSE]

stopifnot(ncol(vst_AG) == 9, ncol(vst_GA) == 9)  # based on your description

# -----------------------------
# 3) Pick soft power in a small-n friendly way
#    With n=9, chasing SFT R^2 > 0.9 often oversparsifies the network.
#    Try a short grid and pick the smallest power giving "acceptable" fit
#    *and* non-collapsed connectivity.
# -----------------------------
pick_power_small_n <- function(datExpr, powers = c(4:12, 14, 16), corFnc = "bicor") {
  # datExpr: samples x genes
  sft <- pickSoftThreshold(
    datExpr,
    powerVector = powers,
    networkType = "signed",
    corFnc = corFnc,
    corOptions = list(maxPOutliers = 0.1),
    verbose = 0
  )
  
  fit <- sft$fitIndices
  # Heuristic: prefer smallest power with SFT.R.sq >= 0.75 and mean.k not tiny
  # (tweak thresholds if needed)
  candidates <- fit[fit$SFT.R.sq >= 0.75 & fit$mean.k. >= 50, , drop = FALSE]
  if (nrow(candidates) == 0) {
    # fallback: maximize (SFT.R.sq) while keeping mean.k reasonable
    # pick the row with highest SFT.R.sq among those with mean.k >= 30, else absolute best SFT
    c2 <- fit[fit$mean.k. >= 30, , drop = FALSE]
    if (nrow(c2) > 0) {
      chosen <- c2[which.max(c2$SFT.R.sq), "Power"]
    } else {
      chosen <- fit[which.max(fit$SFT.R.sq), "Power"]
    }
  } else {
    chosen <- candidates[which.min(candidates$Power), "Power"]
  }
  
  list(power = as.integer(chosen), fitIndices = fit)
}

# MAD per gene within each cross
mad_AG <- apply(vst_AG, 1, mad, constant = 1, na.rm = TRUE)
mad_GA <- apply(vst_GA, 1, mad, constant = 1, na.rm = TRUE)

# also protect against NA
keep <- is.finite(mad_AG) & is.finite(mad_GA) & (mad_AG > 0) & (mad_GA > 0)

vst_shared <- vst_f[keep, ]

# rebuild datExpr matrices (samples x genes)
dat_AG <- t(vst_shared[, coldata$Cross=="AG", drop=FALSE])
dat_GA <- t(vst_shared[, coldata$Cross=="GA", drop=FALSE])

# sanity
stopifnot(all(apply(dat_AG, 2, mad, constant=1, na.rm=TRUE) > 0))
stopifnot(all(apply(dat_GA, 2, mad, constant=1, na.rm=TRUE) > 0))

pAG <- pick_power_small_n(dat_AG)
pGA <- pick_power_small_n(dat_GA)

power_AG <- pAG$power
power_GA <- pGA$power

message("Chosen power AG: ", power_AG)
message("Chosen power GA: ", power_GA)

# Optional: inspect pAG$fitIndices and pGA$fitIndices

# -----------------------------
# 4) Build networks per cross (adjacency + TOM) using bicor, signed
# -----------------------------
make_network <- function(datExpr, power, corFnc = "bicor") {
  # datExpr: samples x genes
  corMat <- cor(
    datExpr,
    use = "pairwise.complete.obs",
    method = if (corFnc == "pearson") "pearson" else "pearson" # bicor handled below
  )
  
  # If using bicor, use WGCNA::bicor for robustness
  if (corFnc == "bicor") {
    corMat <- bicor(datExpr,
                    use = "pairwise.complete.obs",
                    maxPOutliers = 0.1)
  }
  
  adj <- (0.5 * (1 + corMat))^power  # signed adjacency
  diag(adj) <- 0
  
  TOM <- TOMsimilarity(adj, TOMType = "signed")
  dissTOM <- 1 - TOM
  
  list(cor = corMat, adj = adj, TOM = TOM, dissTOM = dissTOM)
}

power <- 9

net_AG <- make_network(dat_AG, power, corFnc = "bicor")
net_GA <- make_network(dat_GA, power, corFnc = "bicor")


# 5) Define modules in ONE reference cross (recommended with n=9)
#    Reference = AG here; you can swap later

geneTree <- hclust(as.dist(net_AG$dissTOM), method = "average")

minModuleSize <- 50
dynamicMods <- cutreeDynamic(
  dendro = geneTree,
  distM = net_AG$dissTOM,
  deepSplit = 2,
  pamRespectsDendro = FALSE,
  minClusterSize = minModuleSize
)

moduleColors <- labels2colors(dynamicMods)

# Merge close modules using eigengenes from reference cross
MEList <- moduleEigengenes(dat_AG, colors = moduleColors, softPower = power_AG)
MEs <- MEList$eigengenes

mergeCutHeight <- 0.25
merge <- mergeCloseModules(dat_AG, moduleColors, cutHeight = mergeCutHeight, verbose = 0)
moduleColors_merged <- merge$colors
MEs_merged <- merge$newMEs

# Final module assignment
ref_colors <- moduleColors_merged
ref_modules <- sort(unique(ref_colors))
ref_modules <- ref_modules[ref_modules != "grey"]  # drop unassigned

message("Reference modules (excluding grey): ", length(ref_modules))

# -----------------------------
# 6) Prepare objects for netrep
#    netrep wants: a list of correlation matrices and a module membership list
# -----------------------------
# Keep only genes not grey
keep <- ref_colors != "grey"
genes_keep <- colnames(dat_AG)[keep]

# Restrict correlation matrices to the same genes
Cref <- net_AG$cor[genes_keep, genes_keep]
Ctest <- net_GA$cor[genes_keep, genes_keep]

# Build module membership list for netrep
modules_list <- split(genes_keep, ref_colors[keep])  # names = module colors
modules_list <- modules_list[names(modules_list) != "grey"]

# Networks as a list (reference then test)
cors <- list(AG = Cref, GA = Ctest)

# -----------------------------
# 7) Run netrep: preservation of AG-defined modules in GA

library(WGCNA)

# Must have identical gene order
stopifnot(identical(colnames(dat_AG), colnames(dat_GA)))
stopifnot(length(ref_colors) == ncol(dat_AG))

# WGCNA "multiData" structure
multiData <- list(
  AG = list(data = dat_AG),
  GA = list(data = dat_GA)
)

# WGCNA "multiColor" structure:
# a list of module assignments; each entry is a list with a color vector per dataset.
# Since modules are defined in AG, you still provide colors for both datasets
# (same vector, since genes are the same and modules refer to genes, not samples).
multiColor <- list(
  AG = ref_colors,
  GA = ref_colors
)



set.seed(1)

pres <- WGCNA::modulePreservation(
  multiData = multiData,
  multiColor = multiColor,
  dataIsExpr = TRUE,
  networkType = "signed",
  corFnc = "bicor",
  corOptions = "maxPOutliers = 0.1, use = 'p'",
  referenceNetworks = 1,        # AG is the reference
  testNetworks = 2,             # GA is the test (optional; can omit and it will test all non-reference)
  nPermutations = 1000, # increase to 1000
  maxModuleSize = 5000,         # OK, but note: your function has maxModuleSize default 1000
  quickCor = 0,
  randomSeed = 1,
  verbose = 3
)

#plot
# See what reference/test labels WGCNA used internally
names(pres$preservation$Z)
# Example outputs: "ref.1", "ref.AG", etc.

# Pick the reference block (you set referenceNetworks=1)
ref_block <- names(pres$preservation$Z)[1]
ref_block
names(pres$preservation$Z[[ref_block]])

test_block <- names(pres$preservation$Z[[ref_block]])[1]
# If multiple exist, choose the one matching your test network by pattern:
test_block <- grep("\\.2$|GA$", names(pres$preservation$Z[[ref_block]]), value = TRUE)[1]

tabZ <- pres$preservation$Z[[ref_block]][[test_block]]
head(tabZ)
Zsummary   <- tabZ[, "Zsummary.pres"]
moduleSize <- tabZ[, "moduleSize"]
module     <- rownames(tabZ)

keep <- module != "grey" & !is.na(Zsummary) & !is.na(moduleSize)

plot(
  moduleSize[keep], Zsummary[keep],
  col = module[keep],
  pch = 19,
  xlab = "Module size",
  ylab = "Zsummary (preservation)",
  main = paste("Module preservation:", ref_block, "->", test_block)
)
abline(h = 2,  lty = 2, col = "blue")
abline(h = 10, lty = 2, col = "darkgreen")
text(moduleSize[keep], Zsummary[keep], labels = module[keep], pos = 3, cex = 0.8)

keep <- rownames(tabZ) != "grey"

plot(tabZ[keep, "Zdensity.pres"], tabZ[keep, "Zconnectivity.pres"],
     col = rownames(tabZ)[keep], pch = 19,
     xlab = "Zdensity.pres",
     ylab = "Zconnectivity.pres",
     main = "Preservation components: density vs connectivity")

abline(h = 0, v = 0, lty = 2, col = "grey50")
text(tabZ[keep, "Zdensity.pres"], tabZ[keep, "Zconnectivity.pres"],
     labels = rownames(tabZ)[keep], pos = 3, cex = 0.8)
# using NetRep
library(NetRep)
# Use the same power you standardized on
power <- 9

# ----- Define the background gene universe for preservation -----
# Option A (recommended): include ALL genes used to build the network/modules, with grey as background
bg_genes <- colnames(dat_AG)  # all genes remaining after your shared filtering

stopifnot(identical(colnames(dat_AG), colnames(dat_GA)))
stopifnot(all(bg_genes %in% colnames(dat_GA)))
stopifnot(length(ref_colors) == length(bg_genes))

# moduleAssignments: named vector (gene -> module label)
moduleAssignments <- ref_colors
names(moduleAssignments) <- bg_genes

# Map grey to background label "0"
moduleAssignments[moduleAssignments == "grey"] <- "0"

# ----- Build correlation list (gene x gene) -----
correlation <- list(
  AG = net_AG$cor[bg_genes, bg_genes],
  GA = net_GA$cor[bg_genes, bg_genes]
)

# ----- Build network list (adjacency; gene x gene) -----
# If you already stored adjacency in net_AG$adj and net_GA$adj, use those.
network <- list(
  AG = net_AG$adj[bg_genes, bg_genes],
  GA = net_GA$adj[bg_genes, bg_genes]
)

# ----- Build data list (samples x genes) -----
data <- list(
  AG = dat_AG[, bg_genes, drop = FALSE],
  GA = dat_GA[, bg_genes, drop = FALSE]
)

# Final sanity checks
stopifnot(identical(colnames(data$AG), colnames(data$GA)))
stopifnot(identical(rownames(correlation$AG), colnames(data$AG)))
stopifnot(identical(colnames(network$AG), colnames(data$AG)))
stopifnot(all(names(moduleAssignments) == colnames(data$AG)))

# ----- Run preservation -----
set.seed(1)
presN <- NetRep::modulePreservation(
  network = network,
  data = data,
  correlation = correlation,
  moduleAssignments = moduleAssignments,
  backgroundLabel = "0",
  discovery = "AG",
  test = "GA",
  nPerm = 1000,
  nThreads = 2,
  null = "overlap",
  alternative = "greater",
  simplify = TRUE,
  verbose = TRUE
)

presN

# explore module AG eigenvalues
#dat_AG is samples x genes (you built it as t(vst_shared[, AG]))
stopifnot(all(rownames(dat_AG) %in% rownames(coldata)))

coldata_AG <- coldata[rownames(dat_AG), , drop = FALSE]

traits_AG <- data.frame(
  GI     = as.numeric(coldata_AG$gi),      # <-- adjust column name if needed (GI vs gi)
  weight = as.numeric(coldata_AG$weight),
  row.names = rownames(coldata_AG)
)

stopifnot(identical(rownames(traits_AG), rownames(dat_AG)))

stopifnot(length(ref_colors) == ncol(dat_AG))
stopifnot(identical(names(ref_colors), colnames(dat_AG)) || is.null(names(ref_colors)))
# If ref_colors has names, ensure they correspond to genes and reorder if needed.

ME_AG <- WGCNA::moduleEigengenes(
  expr   = dat_AG,
  colors = ref_colors
)$eigengenes

ME_AG <- WGCNA::orderMEs(ME_AG)

modTraitCor_AG <- cor(ME_AG, traits_AG, use = "p", method = "pearson")

modTraitP_AG <- WGCNA::corPvalueStudent(
  modTraitCor_AG,
  nSamples = nrow(dat_AG)
)

# FDR across all module×trait tests
modTraitFDR_AG <- matrix(
  p.adjust(modTraitP_AG, method = "BH"),
  nrow = nrow(modTraitP_AG),
  dimnames = dimnames(modTraitP_AG)
)

textMat <- paste0(
  signif(modTraitCor_AG, 2), "\nP=",
  signif(modTraitP_AG, 2)
)
dim(textMat) <- dim(modTraitCor_AG)

WGCNA::labeledHeatmap(
  Matrix = modTraitCor_AG,
  xLabels = colnames(traits_AG),
  yLabels = colnames(ME_AG),
  ySymbols = colnames(ME_AG),
  colorLabels = FALSE,
  colors = WGCNA::blueWhiteRed(50),
  textMatrix = textMat,
  setStdMargins = FALSE,
  cex.text = 0.7,
  zlim = c(-1, 1),
  main = "AG: module–trait relationships"
)

# Example criteria — adjust to taste
sig <- which(abs(modTraitCor_AG) >= 0.3 & modTraitP_AG < 0.05, arr.ind = TRUE)

sig_table <- data.frame(
  module = rownames(modTraitCor_AG)[sig[,1]],
  trait  = colnames(modTraitCor_AG)[sig[,2]],
  cor    = modTraitCor_AG[sig],
  p      = modTraitP_AG[sig],
  FDR    = modTraitFDR_AG[sig]
)

sig_table[order(sig_table$p), ]

# sanity check
mod <- "MEturquoise"  # example
plot(traits_AG$weight, ME_AG[, mod],
     pch = 19, xlab = "weight", ylab = mod,
     main = paste("AG:", mod, "vs weight"))
abline(lm(ME_AG[, mod] ~ traits_AG$weight), col = "red")

# Plotting thinned VCF pure bred ----
## Paths
eigvec_file <- "pbPCA.pca.eigenvec"
eigval_file <- "pbPCA.pca.eigenval"

## Read data
pca <- read.table(eigvec_file, header = FALSE, stringsAsFactors = FALSE)
eigval <- scan(eigval_file)

## Name columns
colnames(pca)[1:2] <- c("FID", "IID")
pc_names <- paste0("PC", seq_len(ncol(pca) - 2))
colnames(pca)[3:ncol(pca)] <- pc_names

pca$group <- ifelse(grepl("^AA", pca$IID), "AA",
                    ifelse(grepl("^GG", pca$IID), "GG", "OTHER"))


var_exp <- eigval / sum(eigval) * 100
pc1_lab <- sprintf("PC1 (%.1f%%)", var_exp[1])
pc2_lab <- sprintf("PC2 (%.1f%%)", var_exp[2])

cols <- c(AA = "firebrick", GG = "steelblue", OTHER = "grey50")

plot(pca$PC1, pca$PC2,
     col = cols[pca$group],
     pch = 19,
     xlab = pc1_lab,
     ylab = pc2_lab,
     main = "PCA of pruned SNPs")

legend("topright",
       legend = names(cols),
       col = cols,
       pch = 19,
       bty = "n")

text(pca$PC1, pca$PC2,
     labels = pca$IID,
     pos = 3,
     cex = 0.7)

plot(pca$PC1, pca$PC3,
     col = cols[pca$group],
     pch = 19,
     xlab = sprintf("PC1 (%.1f%%)", var_exp[1]),
     ylab = sprintf("PC3 (%.1f%%)", var_exp[3]),
     main = "PCA: PC1 vs PC3")

