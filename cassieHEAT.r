#!/usr/bin/env Rscript

### NECESSARY PACKAGES ###

# install CRAN packages
if (!requireNamespace("pheatmap", quietly = TRUE))
  install.packages("pheatmap")
library("pheatmap")

# crea carpeta de resultados si no existe
dir.create("DEG_results", showWarnings = FALSE)

# obtiene una lista de gene comunes
shared_genes <- Reduce(intersect, all_degs)
write.csv(shared_genes, "DEG_results/genes_en_comun.csv", row.names = FALSE)

# empezar una tabla vacía para log2FoldChange
log2fc_matrix <- data.frame(row.names = shared_genes)

# iterar sobre comparaciones para guardar el log2change
for (comp in names(all_degs)) {
  res <- read.csv(paste0("DEG_results/", comp, "_full_results.csv"), row.names = 1)
  # check to see if genes map
  matched <- res[rownames(res) %in% shared_genes, "log2FoldChange"]
  # makes a new column in the datafram
  log2fc_matrix[[comp]] <- matched[match(shared_genes, rownames(res))]
}

# cleans up genes with NA values just in case
log2fc_matrix <- log2fc_matrix[complete.cases(log2fc_matrix), ]

# creates a scaled matrix for better visualization
scaled_matrix <- t(scale(t(log2fc_matrix)))

# Plot the heatmap
pheatmap(scaled_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 6,
         fontsize_col = 10,
         main = "Genes Shared Across Treatments (Log2)")

# saving the heat map
png("DEG_results/heatmap_genes_comunes.png", width = 800, height = 1000, res = 120)
pheatmap(scaled_matrix,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         fontsize_row = 6,
         fontsize_col = 10,
         main = "Genes Shared Across Treatments (Log2)")
dev.off()
