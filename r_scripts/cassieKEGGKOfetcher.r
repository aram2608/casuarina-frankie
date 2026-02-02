# extract unique ko terms from annotation data
unique_kos <- unique(trimws(final_kegg$KEGG_ko))
length(unique_kos)

# map ko terms to kegg pathways with progress and retries
pb <- progress_bar$new(
  format = "[:bar] :current/:total :percent ko: :message",
  total = length(unique_kos), clear = FALSE, width = 60
)

failed_kos <- c()

ko2pathway_list <- lapply(unique_kos, function(ko) {
  pb$tick(tokens = list(message = ko))
  attempt <- 1
  max_retries <- 3
  res <- NULL

  while (attempt <= max_retries) {
    Sys.sleep(0.5)  # respect kegg api rate limits
    res <- tryCatch({
      keggGet(paste0("ko:", ko))[[1]]$PATHWAY
    }, error = function(e) NULL)

    if (!is.null(res)) break
    attempt <- attempt + 1
  }

  if (!is.null(res) && length(res) > 0) {
    data.frame(
      KO = rep(ko, length(res)),
      Pathway_ID = names(res),
      Pathway_Desc = unname(res),
      stringsAsFactors = FALSE
    )
  } else {
    failed_kos <<- c(failed_kos, ko)
    NULL
  }
})

# combine all results into a single dataframe
ko2pathway_df <- do.call(rbind, ko2pathway_list)
head(ko2pathway_df)

# merge ko-pathway mapping with gene annotation
gene2ko <- final_kegg[c("query", "KEGG_ko")]
gene2pathway <- merge(gene2ko, ko2pathway_df, by.x = "KEGG_ko", by.y = "KO")
head(gene2pathway)

# create term2gene and term2name for clusterprofiler
term2gene <- unique(gene2pathway[c("Pathway_ID", "query")])
term2name <- unique(gene2pathway[c("Pathway_ID", "Pathway_Desc")])

# make background gene set
background_kegg <- split(term2gene$query, term2gene$Pathway_ID)