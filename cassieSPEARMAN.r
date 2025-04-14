#!/bin/env Rscript

### Calculates Spearman's Coefft ###

# Extract normalized counts
normalized_counts <- counts(dds, normalized = TRUE)

# Perform Spearman correlation
correlation <- cor.test(normalized_counts[,1], normalized_counts[,2], method = "spearman")

# Print the results
print(correlation)