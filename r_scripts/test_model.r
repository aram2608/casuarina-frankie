suppressPackageStartupMessages(library(DESeq2))
suppressPackageStartupMessages(library(ggplot2))

features_raw <- read.delim("feature_counts.txt", sep = "\t", header = TRUE)

# Keep gene ids + count columns only
cleaned <- features_raw[-c(2, 3, 4, 5, 6)]
counts <- cleaned[, -1]
rownames(counts) <- cleaned[, 1]

coldata <- read.csv("test_coldata.csv", header = TRUE, stringsAsFactors = FALSE)
rownames(coldata) <- coldata$Sample

# Ensure sample order matches count matrix columns
coldata <- coldata[colnames(counts), , drop = FALSE]
stopifnot(identical(rownames(coldata), colnames(counts)))

# Encode predictors as factors
coldata$Condition <- factor(
    coldata$Condition,
    levels = c("control", "cells24", "NINA24", "cells48", "NINA48")
)
coldata$Treatment <- factor(coldata$Treatment, levels = c("control", "cells", "NINA"))
coldata$Time <- factor(coldata$Time, levels = c("0", "24", "48"))
coldata$Replicate <- factor(coldata$Replicate)

# Omnibus test across all observed groups
# LRT is appropriate here because we test "any difference among conditions".
dds_lrt <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = coldata,
    design = ~Condition
)
dds_lrt <- DESeq(dds_lrt, test = "LRT", reduced = ~1, quiet = TRUE)
results_lrt <- results(dds_lrt)

# Time-series interaction test (treated samples only)
# Tests whether time trajectories differ between cells vs NINA.
keep <- coldata$Treatment != "control"
dds_ts <- DESeqDataSetFromMatrix(
    countData = counts[, keep],
    colData = droplevels(coldata[keep, ]),
    design = ~ Treatment + Time + Treatment:Time
)
dds_ts <- DESeq(dds_ts, test = "LRT", reduced = ~ Treatment + Time, quiet = TRUE)
results_ts_interaction <- results(dds_ts)

if (!dir.exists("test_results")) {
    dir.create("test_results", recursive = TRUE)
}

vsd_lrt <- vst(dds_lrt, blind = TRUE)
pca_lrt <- plotPCA(vsd_lrt, intgroup = c("Treatment", "Time"), returnData = TRUE)
var_lrt <- round(100 * attr(pca_lrt, "percentVar"))

plot_lrt <- ggplot(pca_lrt, aes(PC1, PC2, color = Treatment, shape = Time)) +
    geom_point(size = 3, alpha = 0.9) +
    xlab(paste0("PC1: ", var_lrt[1], "% variance")) +
    ylab(paste0("PC2: ", var_lrt[2], "% variance")) +
    ggtitle("PCA - Omnibus LRT model") +
    coord_equal() +
    theme_bw()

ggsave(
    filename = file.path("test_results", "pca_lrt.png"),
    plot = plot_lrt,
    width = 8,
    height = 6,
    dpi = 300
)

vsd_ts <- vst(dds_ts, blind = TRUE)
pca_ts <- plotPCA(vsd_ts, intgroup = c("Treatment", "Time"), returnData = TRUE)
var_ts <- round(100 * attr(pca_ts, "percentVar"))

plot_ts <- ggplot(pca_ts, aes(PC1, PC2, color = Treatment, shape = Time)) +
    geom_point(size = 3, alpha = 0.9) +
    xlab(paste0("PC1: ", var_ts[1], "% variance")) +
    ylab(paste0("PC2: ", var_ts[2], "% variance")) +
    ggtitle("PCA - Time interaction model") +
    coord_equal() +
    theme_bw()

ggsave(
    filename = file.path("test_results", "pca_time_interaction.png"),
    plot = plot_ts,
    width = 8,
    height = 6,
    dpi = 300
)
