#!/usr/bin/env Rscript

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
all_degs <- list()

# loop through pair-wise comparison
for (comp in comparisons) {
  res <- results(dds, name = comp) #saves results for each comparison
  sig <- res[which(res$padj < 0.05 & !is.na(res$padj)), ] # filter for p and NA

  up <- rownames(sig[sig$log2FoldChange > 1.5, ]) # filter for up-reg 1.5
  down <- rownames(sig[sig$log2FoldChange < -1.5, ]) # filter for down-red -1.5

  all_degs[[comp]] <- up # compile up-reg genes
  all_degs[[comp]] <- down # compile down-reg genes
}

# todos los genes arriba y abajo
all_venn <- ggVennDiagram(all_degs, label_alpha = 0, edge_size = 0.5,
                          label = "count",
                          category.names = c("Cells vs. Control 24H",
                                             "Cells vs. Control 48H",
                                             "NINA vs. Control 24H",
                                             "NINA vs. Countrol 48H")) +
  scale_fill_distiller(palette = "Spectral") +
  labs(title = "All Differentially Expressed Genes") +
  theme(text = element_text(size = 14))

# saves the plot??
all_venn

# fixes the long labels
all_venn + scale_x_continuous(expand = expansion(mult = .2))

ggsave("DEG_results/all_DEGs_venn.png",
       plot = all_venn, width = 6, height = 6, dpi = 300)