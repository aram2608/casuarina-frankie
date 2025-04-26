#!/usr/bin/env Rscript

### RUN THIS FIRST AFTER DEG TO ENSURE QUALITY ###

### NECESSARY PACKAGES ###

# install bioconductor packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

# install CRAN packages
if (!requireNamespace("ggplot2", quietly = TRUE))
  install.packages("ggplot2")

# load packages
library("DESeq2")
library("ggplot2")

# creates an output directory
dir.create("DEG_results", showWarnings = FALSE)

# VST transformation
vsd <- vst(dds, blind = TRUE)
pca_data <- plotPCA(vsd, intgroup = c("Treatment"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

# PCA plotting using ggplot2
pca <- ggplot(pca_data, aes(PC1, PC2, color = Treatment)) +
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

# save plot as a png
ggsave("DEG_results/pca_plot.png", plot = pca, width = 6, height = 6, dpi = 300)
