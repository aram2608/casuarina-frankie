if (!require(DESeq2)) {
    suppressPackageStartupMessages(library(DESeq2))
}

# Function to create a deseq object and run DGE analysis
# Returns the deseq object, results() need to be used in order to extract the
# underlying result
deseq <- function(counts, coldata) {
    dds <- DESeqDataSetFromMatrix(
        countData = counts,
        colData = coldata,
        design = ~Condition
    )
    dds <- DESeq(dds, test = "LRT", reduced = ~1, quiet = TRUE)
    return(dds)
}
