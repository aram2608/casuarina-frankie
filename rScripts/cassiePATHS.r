```{r}
# adding annotations for KEGG terms
annotated_C24 <- merge(resC24_df, final_kegg[, c("query", "KEGG_Pathway")], by.x = "gene_id", by.y = "query", all.x = TRUE)
annotated_C48 <-merge(resC48_df, final_kegg[, c("query", "KEGG_Pathway")], by.x = "gene_id", by.y = "query", all.x = TRUE)
annotated_N24 <- merge(resN24_df, final_kegg[, c("query", "KEGG_Pathway")], by.x = "gene_id", by.y = "query", all.x = TRUE)
annotated_N48 <- merge(resN48_df, final_kegg[, c("query", "KEGG_Pathway")], by.x = "gene_id", by.y = "query", all.x = TRUE)

# filtering for genes that have KEGG terms and not "-" values
annotated_C24 <- annotated_C24[!is.na(annotated_C24$KEGG_Pathway) 
                               & annotated_C24$KEGG_Pathway != "-",]

annotated_C48 <- annotated_C48[!is.na(annotated_C48$KEGG_Pathway) 
                               & annotated_C48$KEGG_Pathway != "-",]

annotated_N24 <- annotated_N24[!is.na(annotated_N24$KEGG_Pathway) 
                               & annotated_N24$KEGG_Pathway != "-",]

annotated_N48 <- annotated_N48[!is.na(annotated_N48$KEGG_Pathway) 
                               & annotated_N48$KEGG_Pathway != "-",]

# extract query ids for upregulated genes
upregulated_C24 <- annotated_C24 %>%
  filter(log2FoldChange > 1.5, pvalue < 0.05) %>%
  pull(gene_id)

upregulated_C48 <- annotated_C48 %>%
  filter(log2FoldChange > 1.5, pvalue < 0.05) %>%
  pull(gene_id)

upregulated_N24 <- annotated_N24 %>%
  filter(log2FoldChange > 1.5, pvalue < 0.05) %>%
  pull(gene_id)

upregulated_N48 <- annotated_N48 %>%
  filter(log2FoldChange > 1.5, pvalue < 0.05) %>%
  pull(gene_id)

# now do the same for downregulated genes
downregulated_C24 <- annotated_C24 %>%
  filter(log2FoldChange < -1.5, pvalue < 0.05) %>%
  pull(gene_id)

downregulated_C48 <- annotated_C48 %>%
  filter(log2FoldChange < -1.5, pvalue < 0.05) %>%
  pull(gene_id)

downregulated_N24 <- annotated_N24 %>%
  filter(log2FoldChange < -1.5, pvalue < 0.05) %>%
  pull(gene_id)

downregulated_N48 <- annotated_N48 %>%
  filter(log2FoldChange < -1.5, pvalue < 0.05) %>%
  pull(gene_id)
```


# We should be good to go

```{r}
# enrichment function
enrichment_function <- function(extracted_results, kegg_object) {
  enrichment_table <- enricher(
    extracted_results,
    TERM2GENE = kegg_object,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    minGSSize = 10
  )
  return(enrichment_table)
}

# enrichment function for upregulated
upregulated_e_C24 <- enrichment_function(upregulated_C24, enrich_kegg)
upregulated_e_C24_df <- as.data.frame(upregulated_e_C24)

upregulated_e_C48 <- enrichment_function(upregulated_C48, enrich_kegg)
upregulated_e_C48_df <- as.data.frame(upregulated_e_C48)

upregulated_e_N24 <- enrichment_function(upregulated_N24, enrich_kegg)
upregulated_e_N24_df <- as.data.frame(upregulated_e_N24)

upregulated_e_N48 <- enrichment_function(upregulated_N48, enrich_kegg)
upregulated_e_N48_df <- as.data.frame(upregulated_e_N48)

# enrichment function for downregulated
downregulated_e_C24 <- enrichment_function(downregulated_C24, enrich_kegg)
downregulated_e_C24_df <- as.data.frame(downregulated_e_C24)

downregulated_e_C48 <- enrichment_function(downregulated_C48, enrich_kegg)
downregulated_e_C48_df <- as.data.frame(downregulated_e_C48)

downregulated_e_N24 <- enrichment_function(downregulated_N24, enrich_kegg)
downregulated_e_N24_df <- as.data.frame(downregulated_e_N24)

downregulated_e_N48 <- enrichment_function(downregulated_N48, enrich_kegg)
downregulated_e_N48_df <- as.data.frame(downregulated_e_N48)

# creating an enrichplot function
enrich_obj <- function(results_df, gene_list, background_genes, gene_sets_list) {
  enrichres <- new("enrichResult",
                   readable = FALSE,
                   result = results_df, # dataframe of enrichment results
                   pvalueCutoff = 0.05, # 95% confidence
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05, # 5% of enriched are false positives
                   organism = "UNKNOWN",
                   ontology = "UNKNOWN",
                   gene = gene_list,            # tested gene list
                   keytype = "UNKNOWN",
                   universe = background_genes, # all genes considered
                   gene2Symbol = character(0),
                   geneSets = gene_sets_list)   # mapping from terms to genes
  
  return(enrichres)
}
```

# Creating enrichres objects for plotting using enrichplot

```{r}
# running our functions for each comparson
# c24 up
enrichres_up_C24 <- enrich_obj(results = upregulated_e_C24_df, 
                            gene_list = upregulated_C24, 
                            background_genes = annotation_data$query, # i could extract these and make it
                            gene_sets_list = background_kegg)         # its own dataframe/vector probably

class(enrichres_up_C24)
# c48 up
enrichres_up_C48 <- enrich_obj(results = upregulated_e_C48_df, 
                            gene_list = upregulated_C48, 
                            background_genes = annotation_data$query, 
                            gene_sets_list = background_kegg)

class(enrichres_up_C48)
# n24 up
enrichres_up_N24 <- enrich_obj(results = upregulated_e_N24_df, 
                            gene_list = upregulated_N24, 
                            background_genes = annotation_data$query, 
                            gene_sets_list = background_kegg)

class(enrichres_up_N24)
# n48 up
enrichres_up_N48 <- enrich_obj(results = upregulated_e_N48_df, 
                            gene_list = upregulated_N48, 
                            background_genes = annotation_data$query, 
                            gene_sets_list = background_kegg)

class(enrichres_up_N48)

# plot up
dotplot(enrichres_up_C24) + ggtitle("Upreg Cell vs. Control 24H")
dotplot(enrichres_up_C48) + ggtitle("Upreg Cell vs. Control 48H")
dotplot(enrichres_up_N24) + ggtitle("Upreg NINA vs. Control 24H")
dotplot(enrichres_up_N48) + ggtitle("Upreg NINA vs. Control 48H")
```
```{r}
# running our functions for each comparson
# c24 down
enrichres_down_C24 <- enrich_obj(results = downregulated_e_C24_df, 
                            gene_list = downregulated_C24, 
                            background_genes = annotation_data$query, # i could extract these and make it
                            gene_sets_list = background_kegg)         # its own dataframe/vector probably

class(enrichres_down_C24)
# c48 down
enrichres_down_C48 <- enrich_obj(results = downregulated_e_C48_df, 
                            gene_list = downregulated_C48, 
                            background_genes = annotation_data$query, 
                            gene_sets_list = background_kegg)

class(enrichres_down_C48)
# n24 down
enrichres_down_N24 <- enrich_obj(results = downregulated_e_N24_df, 
                            gene_list = downregulated_N24, 
                            background_genes = annotation_data$query, 
                            gene_sets_list = background_kegg)

class(enrichres_down_N24)
# n48 down
enrichres_down_N48 <- enrich_obj(results = downregulated_e_N48_df, 
                            gene_list = downregulated_N48, 
                            background_genes = annotation_data$query, 
                            gene_sets_list = background_kegg)

class(enrichres_down_N48)
# Plotting down reg pathways
dotplot(enrichres_down_C24) + ggtitle("Downreg Cell vs. Control 24H")
dotplot(enrichres_down_C48) + ggtitle("Downreg Cell vs. Control 48H")
dotplot(enrichres_down_N24) + ggtitle("Downreg NINA vs. Control 24H")
dotplot(enrichres_down_N48) + ggtitle("Downreg NINA vs. Control 48H")
```