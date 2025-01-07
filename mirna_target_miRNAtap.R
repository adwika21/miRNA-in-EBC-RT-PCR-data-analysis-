# Load necessary libraries
library(miRNAtap)

# Define the list of miRNAs
miRNAs <- c('hsa-miR-580*'	,'hsa-miR-92b*')

# Initialize an empty list to store results
all_predictions <- list()

# Loop through each miRNA to predict targets
for (mir in miRNAs) {
  # Get predicted targets for the current miRNA
  predictions <- getPredictedTargets(mir = mir, species = 'hsa', method='geom')
  
  # Convert results to data frame
  predictions_df <- as.data.frame(predictions)
  
  # Add a column for the miRNA ID
  predictions_df$miRNA <- mir
  
  # Append to the list
  all_predictions[[mir]] <- predictions_df
}

# Identify all unique column names across the list
all_colnames <- unique(unlist(lapply(all_predictions, colnames)))

# Standardize columns by filling missing columns with NA
standardize_columns <- function(df, all_colnames) {
  df[setdiff(all_colnames, colnames(df))] <- NA
  df <- df[all_colnames]  # Reorder columns
  return(df)
}

# Apply standardization and combine the standardized data frames
standardized_predictions <- lapply(all_predictions, standardize_columns, all_colnames)
combined_predictions <- do.call(rbind, standardized_predictions)

# Extract the part after the period from target names
clean_target_ID <- sapply(strsplit(rownames(combined_predictions), "\\."), `[`, 2)

# Add cleaned target names as the first column in combined_predictions
updated_combined_predictions <- cbind(target_ID = clean_target_ID, combined_predictions)

# Save the updated data frame to a CSV file
write.csv(updated_combined_predictions, "targets.csv", row.names = FALSE)

