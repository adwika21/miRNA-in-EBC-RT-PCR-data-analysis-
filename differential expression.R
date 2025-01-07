# Load necessary libraries
library(dplyr)
library(tidyr)
library(tidyverse)
library(ggsignif)
library(EnhancedVolcano)

# Read the data
array_data <- read.csv("mirna_copd.csv")

# Prune and clean data
clean_data <- array_data %>%
  mutate_if(is.numeric, ~replace(., . > 35, NA)) %>%
  mutate_if(is.numeric, ~replace(., . <= 18, NA))

# Pivot data to longer format and separate columns
selected_mirna <- clean_data %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Group", values_to = "Ct") %>%
  separate(col = "Group", into = c("Disease_state", "ID")) %>%
  group_by(mature.miRNA.ID) %>%
  summarise(Cnt = n(), NA_len = length(which(is.na(Ct)))) %>%
  ungroup() %>%
  mutate(Non_detect_perc = NA_len / Cnt) %>%
  group_by(mature.miRNA.ID) %>%
  summarise(sel_miRNA = ifelse(Non_detect_perc <= 0.8, 1, 0)) %>%
  ungroup() %>%
  filter(sel_miRNA == 1) %>%
  select(mature.miRNA.ID)

selected_mirna <- as.data.frame(unique(selected_mirna$mature.miRNA.ID))
colnames(selected_mirna) <- "mature.miRNA.ID"

filtered_data <- inner_join(selected_mirna, clean_data, "mature.miRNA.ID") %>%
  select(starts_with(c("mature.miRNA.ID", "Cont", "COPD"))) %>%
  droplevels() %>%
  rename(miRNA_ID = mature.miRNA.ID)

# Normalize data
norm_data <- filtered_data %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Group", values_to = "Ct_prime") %>%
  group_by(Group) %>%
  summarise(Mean_expression = mean(Ct_prime, na.rm = TRUE)) %>%
  ungroup() %>%
  pivot_wider(names_from = Group, values_from = Mean_expression) %>%
  select(c(paste0("Cont.", seq(1, 10, by = 1)), paste0("COPD.", seq(1, 10, by = 1))))

# Convert normalized data to matrix
norm_data <- as.matrix(norm_data[1, ])
norm_data <- matrix(rep(norm_data, each = nrow(data_final)), nrow = nrow(data_final))

# Data normalization
data_final <- data_final %>%
  mutate_if(is.numeric, ~replace_na(., 36)) %>%
  column_to_rownames(var = "miRNA_ID")

data_normalized <- data_final - norm_data
data_normalized <- rownames_to_column(data_normalized, var = "miRNA_ID")

mean_ct_data <- data_normalized %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Group", values_to = "Ct") %>%
  separate(col = "Group", into = c("Disease_state", "ID")) %>%
  mutate(Disease_state = as.factor(Disease_state)) %>%
  group_by(Disease_state, miRNA_ID) %>%
  summarise(Mean_ct = mean(Ct, na.rm = TRUE)) %>%
  ungroup()

fold_change <- mean_ct_data %>%
  pivot_wider(names_from = Disease_state, values_from = Mean_ct) %>%
  group_by(miRNA_ID) %>%
  summarise(fold_change = 2^-(COPD - Cont)) %>%
  ungroup()

selected_mirna <- fold_change %>%
  mutate(diff_exp = case_when(
    fold_change <= 0.5 ~ "DownRegulated",
    fold_change > 0.5 & fold_change < 2 ~ "NoChange",
    fold_change >= 2 ~ "UpRegulated"
  ))

data_input <- data_normalized %>%
  pivot_longer(cols = 2:ncol(.), names_to = "Group", values_to = "Ct") %>%
  separate(col = "Group", into = c("Disease_state", "ID")) %>%
  mutate(Disease_state = as.factor(Disease_state))

ttest_data <- data_input %>%
  group_by(Disease_state, miRNA_ID) %>%
  summarise(value = list(Ct)) %>%
  spread(Disease_state, value) %>%
  group_by(miRNA_ID) %>%
  mutate(p_value = t.test(unlist(Cont), unlist(COPD), paired = FALSE, alternative = "two.sided")$p.value)

ttest_data$padj <- p.adjust(ttest_data$p_value, method = 'fdr')

significant_data <- ttest_data %>%
  inner_join(., selected_mirna, by = "miRNA_ID") %>%
  select(c(miRNA_ID, fold_change, padj, diff_exp)) %>%
  filter(padj < 0.05)

plot_data <- significant_data %>%
  mutate(logFC = log2(fold_change),
         neg_log_padj = -log10(padj),
         significance = case_when(
           logFC > 1 & padj < 0.05 ~ "Upregulated",
           logFC < -1 & padj < 0.05 ~ "Downregulated",
           TRUE ~ "Non-significant"
         ))

# Plot the data
ggplot(plot_data, aes(x = logFC, y = neg_log_padj, color = significance)) +
  geom_point(alpha = 0.5, size = 2) +
  scale_color_manual(values = c('blue', 'red', 'grey')) +
  theme_minimal() +
  geom_vline(xintercept = c(-1, 1), alpha = 0.3) +
  geom_hline(yintercept = -log10(0.05), alpha = 0.3) +
  ggtitle('Control vs COPD') +
  theme(legend.position = "right") +
  labs(color = "Significance", y = "-log10 padj value")

