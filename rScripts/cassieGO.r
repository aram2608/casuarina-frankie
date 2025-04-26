"""
GO Analysis-Prep Part 1

The next bits of code are my attempts at GO analysis. GO terms were found using InterProScan command line tool. The first code chunks are prepping data and creating gene universe.
"""
```{r}
# Import and clean go terms from InterProscan output
raw_goterms <- read.csv("cleaned_gos.tsv", sep = "\t", header = TRUE)

raw_goterms$go_terms <- gsub("InterPro", "", raw_goterms$go_terms)
raw_goterms$go_terms <- gsub("PANTHER", "", raw_goterms$go_terms)
raw_goterms$go_terms <- gsub("\\(\\)", "", raw_goterms$go_terms) # this is some weird regex shenanigans you have to deal with, the \ helps you escape from special characters. So you need one \ to escape from the next backslash, then the next \ helps you escape from the ( and the next two \\ do that same for ). It is really silly but its what chatgpt said to do so whatever
raw_goterms$go_terms <- gsub("\\|", ",", raw_goterms$go_terms)

raw_goterms <- raw_goterms[, c("gene_id", "go_terms")]
write.table(raw_goterms, file = "formatted_gos.tsv", sep = "\t",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
```

"""
# GO Analysis-Prep Part 2

Some of the goterms have redundancy due to PANTHER and Interpro both providing an annotation. To clean this up, I am going to run a python script to only keep the first instance of a goterm for each gene ID.

For example:

gene_ID Go-terms
gene1   GO:1234,GO:1234

Will become just,

gene_ID Go-terms
gene1   GO:1234
"""

```{r}
# Importing the table back in as a gene-2-go object
gene2goMappings <- readMappings(file = "dups_removed_goterms.tsv", sep = "\t", IDsep = ",")
```

# GO Analysis Prep 2

Extracting gene IDs and p-values from data frames with diff expressed genes.

```{r, echo=FALSE,message=FALSE}
# Extract gene IDs and p-values
Cp24 <- dplyr::select(resC24_df, gene_id,padj)
Cp48 <- dplyr::select(resC48_df, gene_id,padj)
Np24 <- dplyr::select(resN24_df, gene_id,padj)
Np48 <- dplyr::select(resN48_df, gene_id,padj)

# sort by p-value
sorted_Cp24 <- Cp24[order(Cp24$padj),]
sorted_Cp48 <- Cp48[order(Cp48$padj),]
sorted_Np24 <- Np24[order(Np24$padj),]
sorted_Np48 <- Np48[order(Np48$padj),]

# convert to topGOs genelist format
topgo_Cp24 <- sorted_Cp24$padj
topgo_Cp48 <- sorted_Cp48$padj
topgo_Np24 <- sorted_Np24$padj
topgo_Np48 <- sorted_Np48$padj

# name the topgo vector
names(topgo_Cp24) <- sorted_Cp24$gene_id
names(topgo_Cp48) <- sorted_Cp48$gene_id
names(topgo_Np24) <- sorted_Np24$gene_id
names(topgo_Np48) <- sorted_Np48$gene_id

# filtering value
padj_cut <- 0.05

# functions for filtering significant genes
top5C24_func <- function(x) {
  return(x < padj_cut)
}

top5C48_func <- function(x) {
  return(x < padj_cut)
}

top5N24_func <- function(x) {
  return(x < padj_cut)
}

top5N48_func <- function(x) {
  return(x < padj_cut)
}
```


# GO analysis

#Building topGOdata object for each Gene Ontology category.

```{r, echo=FALSE, message=FALSE}
# creating GO data objects for molecular function
Cp24_go_mf <- new("topGOdata", 
               ontology = "MF", 
               allGenes = topgo_Cp24, 
               annot = annFUN.gene2GO, 
               geneSel = top5C24_func,
               gene2GO = gene2goMappings)

Cp48_go_mf <- new("topGOdata",
               ontology = "MF",
               allGenes = topgo_Cp48,
               annot = annFUN.gene2GO,
               geneSel = top5C48_func,
               gene2GO = gene2goMappings)

Np24_go_mf <- new("topGOdata",
               ontology = "MF",
               allGenes = topgo_Np24,
               annot = annFUN.gene2GO,
               geneSel = top5N24_func,
               gene2GO = gene2goMappings)

Np48_go_mf <- new("topGOdata",
               ontology = "MF",
               allGenes = topgo_Np48,
               annot = annFUN.gene2GO,
               geneSel = top5N48_func,
               gene2GO = gene2goMappings)

# creating GO data objects for biological process
Cp24_go_bp <- new("topGOdata", 
               ontology = "BP", 
               allGenes = topgo_Cp24, 
               annot = annFUN.gene2GO, 
               geneSel = top5C24_func,
               gene2GO = gene2goMappings)

Cp48_go_bp <- new("topGOdata",
               ontology = "BP",
               allGenes = topgo_Cp48,
               annot = annFUN.gene2GO,
               geneSel = top5C48_func,
               gene2GO = gene2goMappings)

Np24_go_bp <- new("topGOdata",
               ontology = "BP",
               allGenes = topgo_Np24,
               annot = annFUN.gene2GO,
               geneSel = top5N24_func,
               gene2GO = gene2goMappings)

Np48_go_bp <- new("topGOdata",
               ontology = "BP",
               allGenes = topgo_Np48,
               annot = annFUN.gene2GO,
               geneSel = top5N48_func,
               gene2GO = gene2goMappings)

# creating GO data objects for cellular component
Cp24_go_cc <- new("topGOdata", 
               ontology = "CC", 
               allGenes = topgo_Cp24, 
               annot = annFUN.gene2GO, 
               geneSel = top5C24_func,
               gene2GO = gene2goMappings)

Cp48_go_cc <- new("topGOdata",
               ontology = "CC",
               allGenes = topgo_Cp48,
               annot = annFUN.gene2GO,
               geneSel = top5C48_func,
               gene2GO = gene2goMappings)

Np24_go_cc <- new("topGOdata",
               ontology = "CC",
               allGenes = topgo_Np24,
               annot = annFUN.gene2GO,
               geneSel = top5N24_func,
               gene2GO = gene2goMappings)

Np48_go_cc <- new("topGOdata",
               ontology = "CC",
               allGenes = topgo_Np48,
               annot = annFUN.gene2GO,
               geneSel = top5N48_func,
               gene2GO = gene2goMappings)
```
"""
Running Fisher's Test for Enrichment and Plotting

I need to read more into this, the vignette talks about the statistics a bit
but it straight up went over my head. Essentially we are finding which gene functions
are over represented given our gene universe. I also need to find a way to make a function for this, it is so annoying to copy and paste and change values all the time.

Comment from biostars for Ks vs weighted fishers
Fisher and ks are just two ways of answering the same question: are the most significant genes enriched for any particular GO term annotations?

Fisher's exact test compares the expected number of significant genes at random to the observed number of significant genes to arrive at a probability.

The KS test compares the distribution of gene p-values expected at random to the observed distribution of the gene p-values to arrive at a probability. KS is theoretically the better choice because it does not require an arbitrary p-value threshold.

Based of of my most recent project, however, the fisher test with p<0.01 and weight01 algorithm seemed to identify the informative GO terms, whereas KS and weight01 tended to identify very basic GO terms like biological process, or cellular process. This could be particular to my dataset, so I would suggest trying both and seeing which gives you more informative GO terms.
"""

```{r, echo=FALSE,message=FALSE}
# Enrichment tests for Molecular Function Fisher
fisher_Cp24_mf <- runTest(Cp24_go_mf, algorithm = "weight01", statistic = "fisher")
fisher_Cp48_mf <- runTest(Cp48_go_mf, algorithm = "weight01", statistic = "fisher")
fisher_Np24_mf <- runTest(Np24_go_mf, algorithm = "weight01", statistic = "fisher")
fisher_Np48_mf <- runTest(Np48_go_mf, algorithm = "weight01", statistic = "fisher")

# Enrichment tests for Biological Process Fisher
fisher_Cp24_bp <- runTest(Cp24_go_bp, algorithm = "weight01", statistic = "fisher")
fisher_Cp48_bp <- runTest(Cp48_go_bp, algorithm = "weight01", statistic = "fisher")
fisher_Np24_bp <- runTest(Np24_go_bp, algorithm = "weight01", statistic = "fisher")
fisher_Np48_bp <- runTest(Np48_go_bp, algorithm = "weight01", statistic = "fisher")

# Enrichment tests for Cellular Component Fisher
fisher_Cp24_cc <- runTest(Cp24_go_cc, algorithm = "weight01", statistic = "fisher")
fisher_Cp48_cc <- runTest(Cp48_go_cc, algorithm = "weight01", statistic = "fisher")
fisher_Np24_cc <- runTest(Np24_go_cc, algorithm = "weight01", statistic = "fisher")
fisher_Np48_cc <- runTest(Np48_go_cc, algorithm = "weight01", statistic = "fisher")

# Enrichment tests for Molecular Function Kolmogrov-Smirnov
ks_Cp24_mf <- runTest(Cp24_go_mf, algorithm = "classic", statistic = "ks")
ks_Cp48_mf <- runTest(Cp48_go_mf, algorithm = "classic", statistic = "ks")
ks_Np24_mf <- runTest(Np24_go_mf, algorithm = "classic", statistic = "ks")
ks_Np48_mf <- runTest(Np48_go_mf, algorithm = "classic", statistic = "ks")

# Enrichment tests for Biological Process Kolmogrov-Smirnov
ks_Cp24_bp <- runTest(Cp24_go_bp, algorithm = "classic", statistic = "ks")
ks_Cp48_bp <- runTest(Cp48_go_bp, algorithm = "classic", statistic = "ks")
ks_Np24_bp <- runTest(Np24_go_bp, algorithm = "classic", statistic = "ks")
ks_Np48_bp <- runTest(Np48_go_bp, algorithm = "classic", statistic = "ks")

# Enrichment tests for Cellular Component Kolmogrov-Smirnov
ks_Cp24_cc <- runTest(Cp24_go_cc, algorithm = "classic", statistic = "ks")
ks_Cp48_cc <- runTest(Cp48_go_cc, algorithm = "classic", statistic = "ks")
ks_Np24_cc <- runTest(Np24_go_cc, algorithm = "classic", statistic = "ks")
ks_Np48_cc <- runTest(Np48_go_cc, algorithm = "classic", statistic = "ks")

# Enrichment tests for Molecular Function Kolmogrov-Smirnov-Elim
ks_elim_Cp24_mf <- runTest(Cp24_go_mf, algorithm = "elim", statistic = "ks")
ks_elim_Cp48_mf <- runTest(Cp48_go_mf, algorithm = "elim", statistic = "ks")
ks_elim_Np24_mf <- runTest(Np24_go_mf, algorithm = "elim", statistic = "ks")
ks_elim_Np48_mf <- runTest(Np48_go_mf, algorithm = "elim", statistic = "ks")

# Enrichment tests for Biological Process Kolmogrov-Smirnov-Elim
ks_elim_Cp24_bp <- runTest(Cp24_go_bp, algorithm = "elim", statistic = "ks")
ks_elim_Cp48_bp <- runTest(Cp48_go_bp, algorithm = "elim", statistic = "ks")
ks_elim_Np24_bp <- runTest(Np24_go_bp, algorithm = "elim", statistic = "ks")
ks_elim_Np48_bp <- runTest(Np48_go_bp, algorithm = "elim", statistic = "ks")

# Enrichment tests for Cellular Component Kolmogrov-Smirnov-Elim
ks_elim_Cp24_cc <- runTest(Cp24_go_cc, algorithm = "elim", statistic = "ks")
ks_elim_Cp48_cc <- runTest(Cp48_go_cc, algorithm = "elim", statistic = "ks")
ks_elim_Np24_cc <- runTest(Np24_go_cc, algorithm = "elim", statistic = "ks")
ks_elim_Np48_cc <- runTest(Np48_go_cc, algorithm = "elim", statistic = "ks")
```
"""
Generating a Table for Plotting

This function may be broken, need to test it a bit more. I might have to find a different way to do this since the plots look a bit strange.
"""
```{r, echo=FALSE,message=FALSE}
# Function to generate a nicely formatted enrichment table
get_enrich_table <- function(go_obj, fisher_result, ks_result, ks_elim, top_n = 10) {
  tab <- GenTable(go_obj,
                  classicFisher = fisher_result,
                  classicKS = ks_result,
                  elimKS = ks_elim,
                  orderBy = "elimKS",
                  ranksOf = "classicFisher",
                  topNodes = top_n) %>%
    mutate(
      pval = as.numeric(classicFisher),
      elimKS = as.numeric(elimKS),
      classicKS = as.numeric(classicKS),
      Term = factor(Term, levels = rev(Term)),
      Count = as.numeric(Annotated))
  return(tab)
}

# Molecular Function (MF)
tab_Cp24_mf <- get_enrich_table(Cp24_go_mf, 
                                fisher_Cp24_mf, 
                                ks_Cp24_mf, 
                                ks_elim_Cp24_mf)

tab_Cp48_mf <- get_enrich_table(Cp48_go_mf, 
                                fisher_Cp48_mf, 
                                ks_Cp48_mf, 
                                ks_elim_Cp48_mf)

tab_Np24_mf <- get_enrich_table(Np24_go_mf, 
                                fisher_Np24_mf,
                                ks_Np24_mf,
                                ks_elim_Np24_mf)

tab_Np48_mf <- get_enrich_table(Np48_go_mf, 
                                fisher_Np48_mf,
                                ks_Np48_mf,
                                ks_elim_Np48_mf)

# Biological Process (BP)
tab_Cp24_bp <- get_enrich_table(Cp24_go_bp, 
                                fisher_Cp24_bp,
                                ks_Cp24_bp,
                                ks_elim_Cp24_bp)

tab_Cp48_bp <- get_enrich_table(Cp48_go_bp, 
                                fisher_Cp48_bp,
                                ks_Cp48_bp,
                                ks_elim_Cp48_bp)

tab_Np24_bp <- get_enrich_table(Np24_go_bp, 
                                fisher_Np24_bp,
                                ks_Np24_bp,
                                ks_elim_Np24_bp)

tab_Np48_bp <- get_enrich_table(Np48_go_bp, 
                                fisher_Np48_bp,
                                ks_Np48_bp,
                                ks_elim_Np48_bp)

# Cellular Component (CC)
tab_Cp24_cc <- get_enrich_table(Cp24_go_cc, 
                                fisher_Cp24_cc,
                                ks_Cp24_cc,
                                ks_elim_Cp24_cc)

tab_Cp48_cc <- get_enrich_table(Cp48_go_cc, 
                                fisher_Cp48_cc,
                                ks_Cp48_cc,
                                ks_elim_Cp24_cc)

tab_Np24_cc <- get_enrich_table(Np24_go_cc, 
                                fisher_Np24_cc,
                                ks_Np24_cc,
                                ks_elim_Np24_cc)

tab_Np48_cc <- get_enrich_table(Np48_go_cc, 
                                fisher_Np48_cc,
                                ks_Np48_cc,
                                ks_elim_Np24_cc)

"""
# Plotting Bar Plot

Some standard ggplots for enriched GO functions. The p-values are different for each one 
which is why I think something is not quite right in the process. I may need to play around
with the Fisher's test, table, and plot functions to figure out what is going on.
"""
```{r}
# plotting function to include enriched terms as asterisks
plot_enrich_bar <- function(enrich_tab, title) {
  enrich_tab <- enrich_tab %>%
    mutate(asterisk = ifelse(pval < 0.05, "*", ""))

  ggplot(enrich_tab, aes(x = Term, y = Count)) +
    geom_bar(stat = "identity", fill = "#4B9CD3") +
    geom_text(aes(label = asterisk), hjust = -0.5, color = "black", size = 5) +
    coord_flip() +
    labs(title = title,
         x = "GO Term",
         y = "Annotated Genes") +
    theme_minimal()
}

# MF plots
plot_enrich_bar(tab_Cp24_mf, "MF Cells vs. Control 24H")
plot_enrich_bar(tab_Cp48_mf, "MF Cells vs. Control 48H")
plot_enrich_bar(tab_Np24_mf, "MF NINA vs. Control 24H")
plot_enrich_bar(tab_Np48_mf, "MF NINA vs. Control 48H")

# BP plots
plot_enrich_bar(tab_Cp24_bp, "BP Cells vs. Control 24H")
plot_enrich_bar(tab_Cp48_bp, "BP Cells vs. Control 48H")
plot_enrich_bar(tab_Np24_bp, "BP NINA vs. Control 24H")
plot_enrich_bar(tab_Np48_bp, "BP NINA vs. Control 48H")

# CC plots
plot_enrich_bar(tab_Cp24_cc, "CC Cells vs. Control 24H")
plot_enrich_bar(tab_Cp48_cc, "CC Cells vs. Control 48H")
plot_enrich_bar(tab_Np24_cc, "CC NINA vs. Control 24H")
plot_enrich_bar(tab_Np48_cc, "CC NINA vs. Control 48H")
```