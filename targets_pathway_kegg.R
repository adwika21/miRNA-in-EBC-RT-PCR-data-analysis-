#load library
library(ggplot2)
library(org.Hs.eg.db)
library(enrichplot)
library(KEGGREST)
library(clusterProfiler)

target <- read.csv("targets.csv", stringsAsFactors = FALSE)

target <- target$GeneSymbol
head(target)

entrez_ids <- mapIds(org.Hs.eg.db, keys = target, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")

# Remove any NAs that might have resulted from failed mappings
entrez_ids <- entrez_ids[!is.na(entrez_ids)]

# Perform KEGG enrichment with the converted Entrez IDs
kegg_enrich <- enrichKEGG(gene = entrez_ids, organism = 'hsa', pvalueCutoff = 0.05)


# Convert KEGG enrichment results to data frame
kegg_df <- read.csv("kegg_david.csv")

# Plot KEGG enrichment results
ggplot(kegg_df, aes(x = reorder(Term, Count), y = Count)) +
  geom_bar(stat = "identity", fill = "skyblue") +
  coord_flip() +
  theme_minimal() +
  labs(title = "KEGG Pathway Enrichment Analysis", x = "Pathways", y = "Count")

