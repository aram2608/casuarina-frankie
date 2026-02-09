library(ComplexHeatmap)
library(RColorBrewer)
library(jsonlite)

raw_m <- read.delim("matrix.tsv", check.names = FALSE, stringsAsFactors = FALSE)

first_col <- raw_m[[1]]
looks_like_gene_id <- mean(grepl("^g[0-9]+$", first_col)) > 0.8
if (looks_like_gene_id) {
  rownames(raw_m) <- first_col
  raw_m <- raw_m[, -1, drop = FALSE]
}

m <- data.matrix(raw_m)

go_json <- fromJSON("gos_of_interest.json")
genes_of_interest <- names(go_json)

row_group <- rep("Other genes", nrow(m))
if (!is.null(rownames(m))) {
  row_group[rownames(m) %in% genes_of_interest] <- "GO genes"
}

cond <- rep(c("Control", "Frankia 24H", "Supernatant 24H", "Frankia 48H", "Supernatant 48H"), each = 3)

col_labels <- structure(cond, names = colnames(m))

my_palette <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(255)
unique_treatments <- unique(cond)
treatment_palette <- setNames(brewer.pal(length(unique_treatments), "Set1"), unique_treatments)

top_ann <- HeatmapAnnotation(
  Treatment = cond,
  show_annotation_name = FALSE,
  col = list(Treatment = treatment_palette))

right_ann <- rowAnnotation(
  GOI = row_group,
  show_annotation_name = FALSE,
  col = list(GOI = c("GO genes" = "#E64B35", "Other genes" = "#D3D3D3"))
)

map <- Heatmap(m, 
        show_row_names = FALSE, 
        show_column_names = FALSE,
        top_annotation = top_ann,
        right_annotation = right_ann,
        col = my_palette,
        heatmap_legend_param = list(title = "Z-score"))

png(filename = "heatmap.png",
    width = 8,
    height = 6,
    units = "in",
    res = 600)
draw(map)
dev.off()
