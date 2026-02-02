# Use the full unfiltered DESeq2 results (not sorted)
topgo_Cp24 <- resC24_df$padj
topgo_Cp48 <- resC48_df$padj
topgo_Np24 <- resN24_df$padj
topgo_Np48 <- resN48_df$padj

# Set gene names
names(topgo_Cp24) <- resC24_df$gene_id
names(topgo_Cp48) <- resC48_df$gene_id
names(topgo_Np24) <- resN24_df$gene_id
names(topgo_Np48) <- resN48_df$gene_id

# Create named logFC vectors for each comparison
logfc_Cp24 <- resC24_df$log2FoldChange
names(logfc_Cp24) <- resC24_df$gene_id

logfc_Cp48 <- resC48_df$log2FoldChange
names(logfc_Cp48) <- resC48_df$gene_id

logfc_Np24 <- resN24_df$log2FoldChange
names(logfc_Np24) <- resN24_df$gene_id

logfc_Np48 <- resN48_df$log2FoldChange
names(logfc_Np48) <- resN48_df$gene_id

# thresholds for filtering
padj_cutoff <- 0.05
logfc_cutoff <- 1.0

# functions for filtering significant genes
top5C24_func <- function(x) {
  genes <- names(x)
  padj_pass <- x[genes] < padj_cutoff
  fc_pass <- abs(logfc_Cp24[genes]) > logfc_cutoff
  return(padj_pass & fc_pass)
}

top5C48_func <- function(x) {
  genes <- names(x)
  padj_pass <- x[genes] < padj_cutoff
  fc_pass <- abs(logfc_Cp48[genes]) > logfc_cutoff
  return(padj_pass & fc_pass)
}

top5N24_func <- function(x) {
  genes <- names(x)
  padj_pass <- x[genes] < padj_cutoff
  fc_pass <- abs(logfc_Np24[genes]) > logfc_cutoff
  return(padj_pass & fc_pass)
}

top5N48_func <- function(x) {
  genes <- names(x)
  padj_pass <- x[genes] < padj_cutoff
  fc_pass <- abs(logfc_Np48[genes]) > logfc_cutoff
  return(padj_pass & fc_pass)
}