library(multiMiR)
mirnas_lc_cpd <- c("hsa-miR-523", "hsa-miR-922", "hsa-miR-512-3p", "hsa-miR-92b*",
                   "hsa-miR-1911*", "hsa-miR-1185", "hsa-miR-23b*", "hsa-miR-520c-5p",
                   "hsa-miR-95", "hsa-miR-520g", "hsa-miR-296-5p", "hsa-miR-34a",
                   "hsa-miR-527", "hsa-miR-23b")
# Predict targets
miRNA_targets <- get_multimir(mirna = mirnas_lc_cpd, 
                                     table = "validated", 
                                     summary = TRUE)

# Access the data from the mmquery_bioc object
miRNA_targets <- miRNA_targets@data

# Check the first few rows of the data frame
head(miRNA_targets)


# Extract unique miRNA names from the miRNA_targets_df
miRNA_names <- unique(miRNA_targets$mature_mirna_id)

# View the names of the miRNAs
print(miRNA_names)

write.csv(miRNA_targets,"targets.csv")
