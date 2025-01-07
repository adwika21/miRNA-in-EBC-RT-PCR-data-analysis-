# Load necessary libraries
library(dplyr)
library(tidyverse)
library(ggsignif)
library(EnhancedVolcano)

# Read the data
input_data <- read.csv("copd_lc_mirna.csv")

# Clean and process the data
clean_data <- input_data %>%
  mutate_if(is.numeric, ~ifelse(. > 35, NA, .)) %>%
  mutate_if(is.numeric, ~ifelse(. <= 18, NA, .))

# Reshape and filter miRNAs expressed in at least 20% of the samples
miRNA_selection <- clean_data %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Sample", values_to = "Ct_value") %>%
  separate(Sample, into = c("Disease", "SampleID"), sep = "\\.") %>%
  group_by(mature.miRNA.ID) %>%
  summarize(Total_samples = n(), 
            Expressed_samples = sum(!is.na(Ct_value))) %>%
  mutate(Expressed_pct = Expressed_samples / Total_samples) %>%
  filter(Expressed_pct >= 0.5) %>%  # At least 20% of the samples
  select(mature.miRNA.ID)

# Join selected miRNAs back to the data and rename column
filtered_data <- clean_data %>%
  inner_join(miRNA_selection, by = "mature.miRNA.ID") %>%
  select(starts_with(c("mature.miRNA.ID", "copd", "lc"))) %>%
  rename(miRNA = mature.miRNA.ID)

# Normalize data
norm_data <- filtered_data %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Group", values_to = "Ct_prime") %>%
  group_by(Group) %>%
  summarise(Mean_expression = mean(Ct_prime, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(names_from = Group, values_from = Mean_expression) %>%
  select(c(paste0("copd.", seq(1, 40, by = 1)), paste0("lc.", seq(1, 40, by = 1))))

# Convert normalized data to matrix
norm_matrix <- as.matrix(norm_data[1, ])
norm_matrix <- matrix(rep(norm_matrix, each = nrow(filtered_data)), nrow = nrow(filtered_data))

# Data normalization
normalized_data <- filtered_data %>%
  mutate_if(is.numeric, ~replace_na(., 36)) %>%
  column_to_rownames(var = "miRNA") - norm_matrix

normalized_data <- rownames_to_column(normalized_data, var = "miRNA")

# Calculate mean Ct values
mean_ct_data <- normalized_data %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Group", values_to = "Ct") %>%
  separate(col = "Group", into = c("Disease_state", "ID")) %>%
  mutate(Disease_state = as.factor(Disease_state)) %>%
  group_by(Disease_state, miRNA) %>%
  summarise(Mean_ct = mean(Ct, na.rm = TRUE)) %>%
  ungroup()

# Calculate fold change
fold_change <- mean_ct_data %>%
  pivot_wider(names_from = Disease_state, values_from = Mean_ct) %>%
  group_by(miRNA) %>%
  summarise(fold_change = 2^-(lc - copd)) %>%
  ungroup()

# Categorize miRNAs based on differential expression
selected_miRNAs <- fold_change %>%
  mutate(diff_exp = case_when(
    fold_change <= 0.5 ~ "DownRegulated",
    fold_change > 0.5 & fold_change < 2 ~ "NoChange",
    fold_change >= 2 ~ "UpRegulated"
  ))

# Prepare data for t-tests
data_input <- normalized_data %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Group", values_to = "Ct") %>%
  separate(col = "Group", into = c("Disease_state", "ID")) %>%
  mutate(Disease_state = as.factor(Disease_state))

ttest_data <- data_input %>%
  group_by(Disease_state, miRNA) %>%
  summarise(value = list(Ct)) %>%
  spread(Disease_state, value) %>%
  group_by(miRNA) %>%
  mutate(p_value = t.test(unlist(copd), unlist(lc), paired = FALSE, alternative = "two.sided")$p.value)

# Adjust p-values using FDR correction
ttest_data$padj <- p.adjust(ttest_data$p_value, method = 'fdr')

# Identify significant miRNAs
significant_data <- ttest_data %>%
  inner_join(., selected_miRNAs, by = "miRNA") %>%
  select(c(miRNA, fold_change, padj, diff_exp)) %>%
  filter(padj < 0.05)

# Prepare data for plotting
plot_data <- significant_data %>%
  mutate(logFC = log2(fold_change),
         neg_log_padj = -log10(padj),
         significance = case_when(
           padj < 0.05 & logFC > 1 ~ "Upregulated",
           padj < 0.05 & logFC < -1 ~ "Downregulated",
           TRUE ~ "Non-significant"
         ))

# Plot Volcano Plot
ggplot(plot_data, aes(x = logFC, y = neg_log_padj, color = significance)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c('red', 'grey', 'blue')) + 
  theme_minimal() +
  geom_vline(xintercept = c(-1, 1), alpha = 0.3) +  
  geom_hline(yintercept = -log10(0.05), alpha = 0.3) + 
  scale_x_continuous(limits = c(-5, 9), breaks = seq(-5, 9, by = 2)) + 
  scale_y_continuous(limits = c(-1, 9), breaks = seq(-1, 9, by = 2)) +  
  ggtitle('lung cancer vs copd') +  
  theme(legend.position = "right") +
  labs(color = "Significance", y = "-log10 padj")  

# Save the significant data to a CSV file
write.csv(significant_data, file = "significant_miRNA_COPD_vs_TB.csv", row.names = FALSE)
write.csv(plot_data, file = "plot_data_miRNA.csv",row.names = FALSE)

final_output <- plot_data %>%
  select(miRNA,fold_change,padj ,logFC, neg_log_padj, significance)

write.csv(final_output, file = "final_data_miRNA_LC_vs_COPD.csv", row.names = FALSE)

