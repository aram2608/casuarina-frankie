# Pathway Enrichment Analysis using clusterprofiler

Going to be using ClusterProfiler and a tutorial/blog post I found. God bless the bioinformatics community. tidyr used for file cleanup. 

We need to extract genes that have KEGG annotations to make the gene lists and gene universe/background. Otherwise it will water down our analysis which I think happened with the GO term enrichment.

So originally, i just used the raw koterm_example from the $KEGG_Pathway eggnog-annotation, i think it would be better to make my own kegg pathway mappings given the $KO_terms instead, itll give me more control and help catch any pathways that might be missing.

```{r}
# Load in our data
annotation_data <- read.delim("dups_removed_viridae.tsv", sep = "\t", header = TRUE, check.names = FALSE)

# check.names = FALSE helps with R putting X in front of columns that start with a number

# colnames(annotation_data)[1] <- gsub("^X\\.", "", colnames(annotation_data)[1])
# run this if the column names import all funny

# selecting columns of interest
kegg_data <- annotation_data[c(1,13)] # the KEGG_ko column is index 12
kegg_data <- kegg_data[kegg_data$KEGG_Pathway != "-", ]
kegg_data$KEGG_Pathway <- gsub("ko", "map", kegg_data$KEGG_Pathway) # changes ko terms to map terms
final_kegg <- kegg_data
write.csv(final_kegg, "needs_formatting_kegg.tsv", quote = TRUE, row.names = FALSE)
```

# Reformatting

So the pathway enrichment has not been going too well, because KO terms map to a braod array of non relevant pathways, the custom annotations I am trying to make end up enriching for random pathways like cancer prevention, retinol metabolism, etc. Non plant stuff, so I am going to try something different.

In the previous code chunk I reformatted that pathways column KEGG_Pathway so that all the terms start with map. koblank and mapblank are interchangeable when it comes to searching for pathways so either or couldve worked, all I need is consitency accross the terms. Now I am gonna use a python script to collapse duplicates for each column so I don't skew my enrichment and reload the data.

The next part will be trial and error but I am going to gsub out all the non-relevant pathways that got annotated from eggnogg mapper, so human, mouse stuff. That is gonna take some trial and error so we will see how long this takes.

Alright so this is not going to work at alll, I either need to find a better way to annotate KEGG terms or just ignore them for now.

```{r}
# load reformaatted data
final_kegg <- read.delim("formatted_dups_removed_keggs.tsv", sep = "\t", check.names = FALSE, header = TRUE)
```


# Making our custom KO pathway mappings

I had my friend chat help with this, restuctured the old script to take our new aproach. A bunch of boiler plate I did not want to do by hand.

```{r}
# expand KEGG_Pathway field to long format
gene2pathway <- final_kegg %>%
  filter(!is.na(KEGG_Pathway) & KEGG_Pathway != "-") %>%
  select(query, KEGG_Pathway) %>%
  separate_rows(KEGG_Pathway, sep = ",") %>%
  mutate(
    KEGG_Pathway = trimws(KEGG_Pathway),
    Pathway_ID = KEGG_Pathway
  ) %>%
  select(query, Pathway_ID)

# map Pathway_ID → Description (with retries + progress)
unique_maps <- unique(gene2pathway$Pathway_ID)

pb <- progress_bar$new(
  format = "[:bar] :current/:total :percent pathway: :message",
  total = length(unique_maps), clear = FALSE, width = 60
)

failed_maps <- c()
pathway_name_list <- list()

for (map_id in unique_maps) {
  pb$tick(tokens = list(message = map_id))
  attempt <- 1
  max_retries <- 3
  res <- NULL

  while (attempt <= max_retries) {
    Sys.sleep(0.5)
    res <- tryCatch({
      keggFind("pathway", map_id)
    }, error = function(e) NULL)

    if (!is.null(res) && length(res) > 0) break
    attempt <- attempt + 1
  }

  if (!is.null(res) && length(res) > 0) {
    pathway_name_list[[map_id]] <- data.frame(
      Pathway_ID = map_id,
      Pathway_Desc = unname(res),
      stringsAsFactors = FALSE
    )
  } else {
    failed_maps <- c(failed_maps, map_id)
  }
}

term2name <- bind_rows(pathway_name_list)

# final formatting for clusterProfiler
term2gene <- unique(gene2pathway)
term2gene <- term2gene[c(2,1)]

# baackground gene set for custom enrichment
background_kegg <- split(term2gene$query, term2gene$Pathway_ID)

# sanity checks
head(term2gene)
head(term2name)

# alrighty figured it out
# the order for term2gene needs to be Pathway_ID first then query
# so run this after this code chunk
# term2gene <- term2gene(c(2,1))

# then run for a sanity check
# head(term2gene)
```


# Extract genes of interest

Here we need to extract the genes we are interested in, I am picking only the genes that are fold change 1.5 and p-val 0.05. Quite stringent. Dplyr for data wrangling again

```{r, echo=FALSE,message=FALSE}
# trim whitespace for matching
final_kegg$query <- trimws(final_kegg$query)
final_kegg$KEGG_Pathway <- trimws(final_kegg$KEGG_Pathway)
# adding annotations for KEGG terms
kegg_annotated_C24 <- merge(resC24_df, final_kegg[, c("query", "KEGG_Pathway")], by.x = "gene_id", by.y = "query", all.x = TRUE)
kegg_annotated_C48 <-merge(resC48_df, final_kegg[, c("query", "KEGG_Pathway")], by.x = "gene_id", by.y = "query", all.x = TRUE)
kegg_annotated_N24 <- merge(resN24_df, final_kegg[, c("query", "KEGG_Pathway")], by.x = "gene_id", by.y = "query", all.x = TRUE)
kegg_annotated_N48 <- merge(resN48_df, final_kegg[, c("query", "KEGG_Pathway")], by.x = "gene_id", by.y = "query", all.x = TRUE)

# filtering for genes that have KEGG terms and not "-" values
kegg_annotated_C24 <- kegg_annotated_C24[!is.na(kegg_annotated_C24$KEGG_Pathway) 
                               & kegg_annotated_C24$KEGG_Pathway != "-",]

kegg_annotated_C48 <- kegg_annotated_C48[!is.na(kegg_annotated_C48$KEGG_Pathway) 
                               & kegg_annotated_C48$KEGG_Pathway != "-",]

kegg_annotated_N24 <- kegg_annotated_N24[!is.na(kegg_annotated_N24$KEGG_Pathway) 
                               & kegg_annotated_N24$KEGG_Pathway != "-",]

kegg_annotated_N48 <- kegg_annotated_N48[!is.na(kegg_annotated_N48$KEGG_Pathway) 
                               & kegg_annotated_N48$KEGG_Pathway != "-",]

# extract query ids
extracted_C24 <- kegg_annotated_C24 %>%
  filter(abs(log2FoldChange) > 1.5, padj < 0.05) %>%
  pull(gene_id)

extracted_C48 <- kegg_annotated_C48 %>%
  filter(abs(log2FoldChange) > 1.5, padj < 0.05) %>%
  pull(gene_id)

extracted_N24 <- kegg_annotated_N24 %>%
  filter(abs(log2FoldChange) > 1.5, padj < 0.05) %>%
  pull(gene_id)

extracted_N48 <- kegg_annotated_N48 %>%
  filter(abs(log2FoldChange) > 1.5, padj < 0.05) %>%
  pull(gene_id)
```


# We should be good to go

All the mappings are prepped, the GOIs are extracted, it should be good to go for enrichment.

```{r}
# enrichment function
# Enrichment function with TERM2GENE and TERM2NAME as arguments
enrichment_function <- function(extracted_genes, term2gene, term2name) {
  enrichment_table <- enricher(
    gene = extracted_genes,
    TERM2GENE = term2gene,
    TERM2NAME = term2name,
    pvalueCutoff = 0.05,
    pAdjustMethod = "BH",
    qvalueCutoff = 0.05,
    minGSSize = 10
  )
  return(enrichment_table)
}

# running comparisons
enriched_C24 <- enrichment_function(extracted_C24, term2gene, term2name)
enriched_C24_df <- as.data.frame(enriched_C24)

enriched_C48 <- enrichment_function(extracted_C48, term2gene, term2name)
enriched_C48_df <- as.data.frame(enriched_C48)

enriched_N24 <- enrichment_function(extracted_N24, term2gene, term2name)
enriched_N24_df <- as.data.frame(enriched_N24)

enriched_N48 <- enrichment_function(extracted_N48, term2gene, term2name)
enriched_N48_df <- as.data.frame(enriched_N48)

# enrichplot object for plotting
enrich_obj <- function(results_df, gene_list, background_genes, gene_sets_list) {
  enrichres <- new("enrichResult",
                   result = results_df,
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   qvalueCutoff = 0.05,
                   gene = gene_list,
                   universe = background_genes,
                   geneSets = gene_sets_list,
                   keytype = "UNKNOWN",
                   organism = "ko",
                   ontology = "KEGG",
                   readable = FALSE,
                   gene2Symbol = setNames(gene_list, gene_list)) # fallback mapping
  return(enrichres)
}
```

# Creating enrichres objects for plotting using enrichplot

```{r}
# running our functions for each comparson
enrichres_C24 <- enrich_obj(results_df = enriched_C24_df, 
                            gene_list = extracted_C24, 
                            background_genes = annotation_data$query, # all the genes in our RNA seq experiment
                            gene_sets_list = background_kegg) # custom gene set we made earlier

class(enrichres_C24)
# c48
enrichres_C48 <- enrich_obj(results_df = enriched_C48_df, 
                            gene_list = extracted_C48, 
                            background_genes = annotation_data$query, 
                            gene_sets_list = background_kegg)

class(enrichres_C48)
# n24
enrichres_N24 <- enrich_obj(results_df = enriched_N24_df, 
                            gene_list = extracted_N24, 
                            background_genes = annotation_data$query, 
                            gene_sets_list = background_kegg)

class(enrichres_N24)
# n48
enrichres_N48 <- enrich_obj(results_df = enriched_N48_df, 
                            gene_list = extracted_N48, 
                            background_genes = annotation_data$query, 
                            gene_sets_list = background_kegg)

class(enrichres_N48)

# plot
C24_dotplot <- dotplot(enrichres_C24) + ggtitle("Cell vs. Control 24H")
C48_dotplot <- dotplot(enrichres_C48) + ggtitle("Cell vs. Control 48H")
N24_dotplot <- dotplot(enrichres_N24) + ggtitle("Supernatant vs. Control 24H")
N48_dotplot <- dotplot(enrichres_N48) + ggtitle("Supernatant vs. Control 48H")

# saving plots
ggsave("DEG_results/dotplot_kegg_C24.png",
       plot = C24_dotplot, width = 8, height = 6, dpi = 300)
ggsave("DEG_results/dotplot_kegg_C48.png",
       plot = C48_dotplot, width = 8, height = 6, dpi = 300)
ggsave("DEG_results/dotplot_kegg_N24.png",
       plot = N24_dotplot, width = 8, height = 6, dpi = 300)
ggsave("DEG_results/dotplot_kegg_N48.png",
       plot = N48_dotplot, width = 8, height = 6, dpi = 300)

# plot to markdown
C24_dotplot
C48_dotplot
N24_dotplot
N48_dotplot
```