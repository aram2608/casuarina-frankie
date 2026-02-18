if (!require(DESeq2)) {
    suppressPackageStartupMessages(library(DESeq2))
}

if (!require(ggplot2)) {
    suppressPackageStartupMessages(library(ggplot2))
}

pca_plot_lrt <- function(dds_obj) {
    vsd_lrt <- vst(dds_obj, blind = TRUE)
    pca_lrt <- plotPCA(vsd_lrt, intgroup = c("Treatment", "Time"), returnData = TRUE)
    var_lrt <- round(100 * attr(pca_lrt, "percentVar"))

    plot_lrt <- ggplot(pca_lrt, aes(PC1, PC2, color = Treatment, shape = Time)) +
        geom_point(size = 3, alpha = 0.9) +
        xlab(paste0("PC1: ", var_lrt[1], "% variance")) +
        ylab(paste0("PC2: ", var_lrt[2], "% variance")) +
        coord_fixed() +
        ggtitle("Omnibus LRT model") +
      scale_color_brewer(palette = "Set1") +
      theme_minimal(base_size = 14) +
      theme(
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(color = "black"),
        panel.grid.major = element_line(color = "gray90"),
        panel.grid.minor = element_blank()
      )
    return(plot_lrt)
}
