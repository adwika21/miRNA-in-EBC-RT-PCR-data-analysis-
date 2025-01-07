library(ReactomePA)

reactome_enrich <- enrichPathway(gene = entrez_ids, organism = 'human', pvalueCutoff = 0.05)

# Convert Reactome enrichment results to data frame
reactome_df <- as.data.frame(reactome_enrich)

# Plot Reactome enrichment results
ggplot(reactome_df, aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", fill = "lightcoral") +
  coord_flip() +
  theme_minimal() +
  labs(title = "Reactome Pathway Enrichment Analysis", x = "Pathways", y = "-log10(p-value)")

