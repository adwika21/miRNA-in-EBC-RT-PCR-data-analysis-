
library(org.Hs.eg.db)
# Load your target genes
target_genes <- read.csv("targets.csv", stringsAsFactors = FALSE)


# Convert Entrez IDs to Gene Symbols
entrez_ids <- target_genes$target_ID
# Convert entrez_ids to character vector
entrez_ids <- as.character(target_genes$target_ID)

# Now map the Entrez IDs to gene symbols
gene_symbols <- mapIds(org.Hs.eg.db, keys = entrez_ids, column = "SYMBOL", keytype = "ENTREZID", multiVals = "first")


# Add the gene symbols to your target data frame
target_genes$GeneSymbol <- gene_symbols

# Remove rows where gene symbol could not be mapped (NA values)
target_genes <- target_genes[!is.na(target_genes$GeneSymbol), ]

# Save the updated target data with gene symbols
write.csv(target_genes, file = "targets_symbol.csv", row.names = FALSE)

