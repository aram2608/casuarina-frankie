#!/usr/bin/env Rscript

# Function to load in counts file, performs a perfunctory cleaning
# step and sets rownames
load_counts <- function(file_path) {
    raw <- read.delim(file_path, sep = "\t", header = TRUE)
    cleaned <- raw[-c(2:6)]
    counts <- cleaned[, -1]
    rownames(counts) <- cleaned[, 1]
    colnames(counts) <- trimws(colnames(counts))
    return(counts)
}

# Loads in the column data
load_coldata <- function(file_path) {
  coldata <- read.csv(file_path, header = TRUE, stringsAsFactors = FALSE)
  coldata$Sample <- trimws(coldata$Sample)
  rownames(coldata) <- coldata$Sample
  return(coldata)
}
