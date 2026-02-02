#!/usr/bin/env Rscript

### RUN THIS AFTER PCA ###
### CREATES NECESSARY DATA FRAMES FOR FURTHER PLOTS ###

### NECESSARY PACKAGES ###

# installing CRAN packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("ggVennDiagram", quietly = TRUE))
  install.packages("ggVennDiagram")
if (!requireNamespace("ggplo2", quietly = TRUE))
  install.packages("ggplot2")
# install bioconductor packages
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

# load libraries
library("DESeq2")
library("ggVennDiagram")
library("ggplot2")

# creates an output directory
dir.create("DEG_results", showWarnings = FALSE)

# list to store genes
comparisons <- resultsNames(dds)[-1] # removes the first weird column
upregulated_genes <- list()
downregulated_genes <- list()
all_degs <- list()

# loop through pair-wise comparison
for (comp in comparisons) {
  res <- results(dds, name = comp) #saves results for each comparison
  sig <- res[which(res$padj < 0.05 & !is.na(res$padj)), ] # filter for p-val

  up <- rownames(sig[sig$log2FoldChange > 1.5, ]) # filter for up-reg 1.5
  down <- rownames(sig[sig$log2FoldChange < -1.5, ]) # filter for down-red -1.5
  all <- rownames(sig) # collection of DEGs

  upregulated_genes[[comp]] <- up # compile upregulated
  downregulated_genes[[comp]] <- down # compile downregulated
  all_degs[[comp]] <- all # compile all DEGs
}

# plot diagrams for up, down, and all genes
# upreg genes
ggVennDiagram(upregulated_genes, label_alpha = 0, edge_size = 0.5) +
  scale_fill_gradient(low = "yellow", high = "red") +
  labs(title = "Upregulated Genes") +
  theme(text = element_text(size = 14))
ggsave("DEG_results/upregulated_venn.png", width = 6, height = 6, dpi = 300)

# downreg genes
ggVennDiagram(downregulated_genes, label_alpha = 0, edge_size = 0.5) +
  scale_fill_gradient(low = "purple", high = "blue") +
  labs(title = "Downregulated Genes") +
  theme(text = element_text(size = 14))
ggsave("DEG_results/downregulated_venn.png", width = 6, height = 6, dpi = 300)

# todos los genes arriba y abajo
ggVennDiagram(all_degs, label_alpha = 0, edge_size = 0.5) +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  labs(title = "All Differentially Expressed Genes") +
  theme(text = element_text(size = 14))
ggsave("DEG_results/all_DEGs_venn.png", width = 6, height = 6, dpi = 300)
