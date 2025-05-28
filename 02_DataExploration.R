rm(list = ls())

# Load libraries
library(tidyverse)
library(DESeq2)
library(Hmisc)
library(pheatmap)
library(vsn)

# Folder paths
tidy_data_dir <- "01_TidyData"
qc_dir <- "02_QualityCheck"
if (!dir.exists(qc_dir)) dir.create(qc_dir)

# Read data
all_labels <- read.csv(file.path(tidy_data_dir, "all_strains_labels.csv")) %>%
  mutate(condition = paste(cSource, time, sep = "_"))

all_counts <- read.csv(file.path(tidy_data_dir, "all_strains_counts.csv"))

# Reference levels (editable)
ref_cSource <- "glucose"         # <-- your actual reference cSource name
ref_time <- "5d"                 # <-- your actual reference time name
ref_condition <- "glucose_5d"    # <-- your actual reference condition name
ref_strain <- "CBS_464"          # <-- your actual reference strain name

# Prepare count matrix
counts_matrix <- all_counts %>%
  column_to_rownames("protein_id") %>%
  as.matrix()

all_labels <- all_labels %>%
  filter(sample %in% colnames(counts_matrix)) %>%
  column_to_rownames("sample")

counts_matrix <- counts_matrix[, rownames(all_labels)]

# Factor conversion and releveling
all_labels <- all_labels %>%
  mutate(
    condition = factor(condition),
    strain = factor(strain),
    time = factor(time),
    cSource = factor(cSource)
  )

all_labels$cSource   <- relevel(all_labels$cSource, ref = ref_cSource)
all_labels$time      <- relevel(all_labels$time, ref = ref_time)
all_labels$condition <- relevel(all_labels$condition, ref = ref_condition)
all_labels$strain    <- relevel(all_labels$strain, ref = ref_strain)

# Design formula (editable)
# Example: design_formula <- ~ strain + cSource + time + condition
design_formula <- ~ strain + condition

# DESeq2 pipeline
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_matrix),
  colData = all_labels,
  design = design_formula
)

dds <- DESeq(dds)
norm_counts <- counts(dds, normalized = TRUE)

# Pearson correlation
cor_result <- rcorr(as.matrix(norm_counts), type = "pearson")
cor_matrix <- cor_result$r
p_matrix   <- cor_result$P

# Flag outlier replicates
threshold <- 0.8

meta_with_ids <- all_labels %>%
  rownames_to_column("sample") %>%
  mutate(group_id = paste(cSource, time, strain, condition, sep = "_"))

# Long format correlation table
cor_long <- as.data.frame(cor_matrix) %>%
  rownames_to_column("sample1") %>%
  pivot_longer(-sample1, names_to = "sample2", values_to = "correlation") %>%
  filter(sample1 != sample2)  # Exclude self-correlations

# Add group IDs
cor_long <- cor_long %>%
  left_join(meta_with_ids %>% select(sample1 = sample, group_id1 = group_id), by = "sample1") %>%
  left_join(meta_with_ids %>% select(sample2 = sample, group_id2 = group_id), by = "sample2") %>%
  filter(group_id1 == group_id2) %>%
  select(group_id = group_id1, sample1, sample2, correlation)

# Flag groups with ANY low correlation
flagged_groups <- cor_long %>%
  filter(correlation < threshold) %>%
  distinct(group_id)

# Mark replicates with low correlation to most others
outlier_samples <- cor_long %>%
  filter(group_id %in% flagged_groups$group_id) %>%
  group_by(group_id, sample1) %>%
  summarise(
    total_peers = n(),
    n_low_corr = sum(correlation < threshold),
    .groups = "drop"
  ) %>%
  filter(n_low_corr > 1)  # More than 1 low correlation is safer than requiring *all*

# Output results
print("Groups with low correlation pairs:")
print(flagged_groups)

print("Samples with low correlation to all peers:")
print(outlier_samples)

# Save results
write_csv(flagged_groups, file.path(qc_dir, "flagged_replicate_groups.csv"))
write_csv(outlier_samples, file.path(qc_dir, "replicates_low_corr_to_all_peers.csv"))
write.csv(cor_matrix, file.path(qc_dir, "pearson_correlation_matrix.csv"))
write.csv(p_matrix, file.path(qc_dir, "pearson_pvalue_matrix.csv"))

# Save heatmap
pheatmap(cor_matrix,
         main = "Sample Pearson Correlation (Normalized Counts)",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         display_numbers = FALSE,
         filename = file.path(qc_dir, "pearson_correlation_heatmap.pdf"))

# Remove outlier samples from metadata and counts
if (nrow(outlier_samples) > 0) {
  outlier_ids <- outlier_samples$sample1
  all_labels_clean <- all_labels[!rownames(all_labels) %in% outlier_ids, ]
  counts_matrix_clean <- counts_matrix[, !colnames(counts_matrix) %in% outlier_ids]
} else {
  message("No outliers found. Using full dataset.")
  all_labels_clean <- all_labels
  counts_matrix_clean <- counts_matrix
}

# Re-create DESeq2 dataset with cleaned data
dds_clean <- DESeqDataSetFromMatrix(
  countData = round(counts_matrix_clean),
  colData = all_labels_clean,
  design = design_formula
)

# Run DESeq2 on cleaned data
dds_clean <- DESeq(dds_clean)

# Calculate the different transformation methods (blind since we are texting data quality)
ntd <- normTransform(dds_clean)
vsd <- vst(dds_clean, blind = T)
rld <- rlog(dds_clean, blind = T)

# look at the mean SD plots to select a transformation method
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))
meanSdPlot(assay(rld))

# Select the smoothest transformation
dds_transformed <- rld

# Filter the genes based on counts > 10
selected_genes <- ntd %>% assay %>% data.frame() %>%
  filter(rowSums(across(everything(), ~ . > 10)) > 0) %>%
  rownames_to_column(var = "protein_id")

# Create the dds and transformed with the corresponding genes
selectedgenecols <- which(rownames(dds_clean) %in% selected_genes$protein_id) %>%
  unique()

dds_clean <- dds[selectedgenecols, ]
dds_transformed <- dds_transformed[selectedgenecols, ]

# Perform PCA (edit to display the groupings you want)
pca_data <- plotPCA(dds_transformed, intgroup = c("strain", "condition"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

# Save PCA plot
pca_plot <- ggplot(pca_data, aes(PC1, PC2, color = condition, shape = strain)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal(base_size = 14) +
  ggtitle("PCA of Samples (Outliers Removed)")

ggsave(filename = file.path(qc_dir, "pca_plot_cleaned.pdf"), plot = pca_plot, width = 8, height = 6)

# Optional: print plot to viewer
print(pca_plot)

