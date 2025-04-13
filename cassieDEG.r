#!/usr/bin/env Rscript

### THESE VALUES ARE HARD CODED ###
### MAKE SURE TO CHANGE IF RE-USING ###

### REQUIRED PACKAGES ###

# install bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if (!requireNamespace("DESEq2", quietly = TRUE))
  BiocManager::install("DESeq2")

# load libraries
library("DESeq2")

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