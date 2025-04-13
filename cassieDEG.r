#!/usr/bin/env Rscript

### THESE VALUES ARE HARD CODED ###
### MAKE SURE TO CHANGE IF RE-USING ###

### REQUIRED PACKAGES ###
# install bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
# install CRAN packages
install.packages("ggplot2")
install.packages("pheatmap")
# load libraries
library("DESeq2")
library("ggplot2") #ggplot2 tidyr tibble
library("pheatmap")

### IMPORTING DATA ###
# analyzing count matrix
featurecounts_raw <- read.delim("feature_counts.txt",
                                sep = "\t",
                                header = TRUE) # import feature_counts
cleanfeatures <- featurecounts_raw[-c(2, 3, 4, 5, 6)] # clean up columns
# set column names
counts <- cleanfeatures[, -1]
rownames(counts) <- cleanfeatures[, 1]
# import column data
coldata <- read.csv("coldata.csv", header = TRUE)
# set row names for coldata
rownames(coldata) <- coldata$X
# test to see if columns and rows match
all(rownames(coldata) == colnames(counts)) # should be true

### DEG ANALYSIS ###
# creates a DeseqDataSet or DDS
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ Treatment)
# sets the comparisons, in this case the coldata column name is Treatment,
# the reference is named control
dds$Treatment <- relevel(dds$Treatment, ref = "control")
# no clue what this does, manual said to tho
# prints the dataframe i suppose?
dds
# performs DEG
dds <- DESeq(dds)
# creates a results dataframe
results <- results(dds)
# prints results
results
# prints summary of results
summary(results)
# prints pairwise comparisons,useful double check
resultsNames(dds)

# PCA PLOT
# transformation using VST
vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = c("Treatment"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

# plotting using ggplot
ggplot(pca_data, aes(PC1, PC2, color = Treatment)) +
  geom_point(size = 3, alpha = 0.8) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  coord_fixed() +
  scale_color_brewer(palette = "Set1") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.major = element_line(color = "gray90"),
    panel.grid.minor = element_blank()
  )

# plots an MA plot, another way to see quality
plotMA(results, ylim = c(-2, 2))
