#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(argparse)
  library(KEGGREST)
  library(tidyr)
  library(dplyr)
})

# -------------------------------
# Argument parser
# -------------------------------
parser <- ArgumentParser(description = "KO → KEGG Pathway Mapper")

parser$add_argument("-i", "--input", required = TRUE,
                    help = "Input annotation file (e.g. glauca_annotations_clean.txt)")
parser$add_argument("-c", "--column", required = TRUE, type = "integer",
                    help = "Column index (1-based) for KO terms (e.g. 12)")
parser$add_argument("-o", "--output", default = "ko2pathway_df.rds",
                    help = "Output RDS file for KO to pathway mapping")

args <- parser$parse_args()

input_file <- args$input
column_index <- args$column
output_file <- args$output

cat("Reading:", input_file, "\n")
annotation_data <- read.delim(input_file, sep = "\t", header = TRUE)

# -------------------------------
# Extract and clean KO terms
# -------------------------------
cat("Extracting KO terms from column index:", column_index, "\n")
kegg_data <- annotation_data[, c(1, column_index)]
colnames(kegg_data) <- c("query", "KEGG_ko")

kegg_data <- kegg_data[kegg_data$KEGG_ko != "-", ]
kegg_data$KEGG_ko <- gsub("ko:", "", kegg_data$KEGG_ko)

# Keep all KO terms, not just the first
final_kegg <- separate_rows(kegg_data, KEGG_ko, sep = ",")

# -------------------------------
# Unique KO terms
# -------------------------------
unique_kos <- unique(trimws(final_kegg$KEGG_ko))
unique_kos <- unique_kos[grepl("^K\\d{5}$", unique_kos)]

cat("Unique KO terms to query:", length(unique_kos), "\n")

# -------------------------------
# Load previous session if exists
# -------------------------------
if (file.exists(output_file)) {
  cat("Resuming from saved file...\n")
  ko2pathway_list <- readRDS(output_file)
} else {
  ko2pathway_list <- list()
}

remaining_kos <- setdiff(unique_kos, names(ko2pathway_list))

# -------------------------------
# Query KEGG
# -------------------------------
for (i in seq_along(remaining_kos)) {
  ko <- remaining_kos[i]
  cat("[", i, "/", length(remaining_kos), "] Querying:", ko, "\n")

  res <- tryCatch({
    keggGet(paste0("ko:", ko))[[1]]$PATHWAY
  }, error = function(e) NULL)

  if (!is.null(res)) {
    ko2pathway_list[[ko]] <- data.frame(
      KO = rep(ko, length(res)),
      Pathway_ID = names(res),
      Pathway_Desc = unname(res),
      stringsAsFactors = FALSE
    )
  } else {
    ko2pathway_list[[ko]] <- NULL
  }

  if (i %% 10 == 0 || i == length(remaining_kos)) {
    saveRDS(ko2pathway_list, output_file)
  }

  Sys.sleep(0.5)
}

cat("KEGG mapping complete. Writing combined data frame to:", output_file, "\n")
ko2pathway_df <- do.call(rbind, ko2pathway_list)
saveRDS(ko2pathway_df, output_file)
cat("Done. Total mapped KO terms:", nrow(ko2pathway_df), "\n")

# example usage
"
chmod +x kegg_mapper.R  # make executable (optional)
Rscript kegg_mapper.R \
  -i glauca_annotations_clean.txt \
  -c 12 \
  -o ko2pathway_df.rds

"