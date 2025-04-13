#!/usr/bin/env Rscript

### RUN THIS AFTER PCA ###
### CREATES NECESSARY DATA FRAMES FOR FURTHER PLOTS ###

### PAQUETES NECESARIOS ###

# installing CRAN packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("ggVennDiagram", quietly = TRUE))
  install.packages("ggVennDiagram")
# install bioconductor packages
if (!requireNamespace("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")

# load libraries
library("DESeq2")
library("ggVennDiagram")

# creates an output directory
dir.create("DEG_results", showWarnings = FALSE)

# list to store genes
comparisons <- resultsNames(dds)[-1]

upregulated_genes <- list()
downregulated_genes <- list()
all_degs <- list()

# loop through comparison
for (comp in comparisons) {
  res <- results(dds, name = comp)
  sig <- res[which(res$padj < 0.05 & !is.na(res$padj)), ]

  up <- rownames(sig[sig$log2FoldChange > 0, ])
  down <- rownames(sig[sig$log2FoldChange < 0, ])
  all <- rownames(sig)

  upregulated_genes[[comp]] <- up
  downregulated_genes[[comp]] <- down
  all_degs[[comp]] <- all

  # Exportar resultados
  write.csv(as.data.frame(res), file = paste0("DEG_results/", comp, "_full_results.csv"))
  write.csv(up, file = paste0("DEG_results/", comp, "_upregulated.csv"), row.names = FALSE)
  write.csv(down, file = paste0("DEG_results/", comp, "_downregulated.csv"), row.names = FALSE)
  write.csv(all, file = paste0("DEG_results/", comp, "_all_sig.csv"), row.names = FALSE)
}

# plot diagrams for up, down, and all genes

# upreg genes
ggVennDiagram(upregulated_genes, label_alpha = 0, edge_size = 0.5) +
  scale_fill_gradient(low = "white", high = "red") +
  labs(title = "Genes Upregulated (padj < 0.05, log2FC > 0)") +
  theme(text = element_text(size = 14))
ggsave("DEG_results/upregulated_venn.png", width = 6, height = 6, dpi = 300)

# downreg genes
ggVennDiagram(downregulated_genes, label_alpha = 0, edge_size = 0.5) +
  scale_fill_gradient(low = "white", high = "blue") +
  labs(title = "Genes Downregulated (padj < 0.05, log2FC < 0)") +
  theme(text = element_text(size = 14))
ggsave("DEG_results/downregulated_venn.png", width = 6, height = 6, dpi = 300)

# todos los genes arriba y abajo
ggVennDiagram(all_degs, label_alpha = 0, edge_size = 0.5) +
  scale_fill_gradient(low = "white", high = "darkgreen") +
  labs(title = "Todos los genes diferencialmente expresados (padj < 0.05)") +
  theme(text = element_text(size = 14))
ggsave("DEG_results/all_DEGs_venn.png", width = 6, height = 6, dpi = 300)
