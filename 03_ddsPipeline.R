rm(list = ls())

# Load libraries
library(tidyverse)
library(DESeq2)

# Folder paths
tidy_data_dir <- "01_TidyData"
qc_dir <- "02_QualityCheck"
dds_dir <- "03_DESeqObjects"
if (!dir.exists(dds_dir)) dir.create(dds_dir)

# Load data
all_labels <- read.csv(file.path(tidy_data_dir, "all_strains_labels.csv")) %>%
  mutate(condition = paste(cSource, time, sep = "_"),
         strain_condition = paste(strain, condition, sep = "_"))

all_counts <- read.csv(file.path(tidy_data_dir, "all_strains_counts.csv"))

# Load outliers (from previous script)
outlier_file <- file.path(qc_dir, "replicates_low_corr_to_all_peers.csv")

# If outliers exist, read and remove from data
if (file.exists(outlier_file)) {
  outliers <- read.csv(outlier_file)
  outlier_ids <- outliers$sample1
  message("Outliers detected and will be removed:")
  print(outlier_ids)
  
  all_labels <- all_labels %>%
    filter(!sample %in% outlier_ids)
  
  all_counts <- all_counts %>%
    select(-all_of(outlier_ids))
} else {
  message("No outlier file found — using full dataset.")
}

# Define the strains and conditions you want to keep
strains_to_keep <- unique(all_labels$strain)       # Update with your actual strain names
conditions_to_keep <- unique(all_labels$condition) # Update with the relevant conditions

# Filter the metadata and counts to keep
all_labels <- all_labels %>%
  filter(strain %in% strains_to_keep,
         condition %in% conditions_to_keep)

all_counts <- all_counts %>%
  select(protein_id, all_of(all_labels$sample))

# Define strain-condition reference pairs
strain_condition_combos <- all_labels %>%
  select(strain, condition, strain_condition) %>%
  distinct()

# Tracker for rlog, VST and norm saves
vst_norm_saved <- list()

#i = 26 # for debugging

# Run pipeline for each (strain_ref, cond_ref) pair
for (i in 1:nrow(strain_condition_combos)) {
  strain_ref <- strain_condition_combos$strain[i]
  cond_ref   <- strain_condition_combos$condition[i]
  strain_cond_ref   <- strain_condition_combos$strain_condition[i]
  
  #️ Status
  message(paste0("\n===== Reference: ", strain_ref, " / ", cond_ref, " ====="))
  
  # Conditions present in this strain only
  strain_ref_conditions <- strain_condition_combos %>%
    filter(strain == strain_ref) %>%
    pull(condition)
  
  # Define the strain sets to analyze
  strain_sets <- list(
    single_strain = list(strain_ref),
    all_strains = as.list(unique(all_labels$strain))
  )
  
  #strain_set_name <- "all_strains" # for debugging
  
  for (strain_set_name in names(strain_sets)) {
    strain_group <- unlist(strain_sets[[strain_set_name]])
    
    message(paste0("→ Processing: ", strain_set_name, " (strain = ", strain_ref, ", condition = ", cond_ref, ")"))
    
    # Subset metadata and counts
    dds_labels <- all_labels %>%
      filter(strain %in% strain_group, condition %in% strain_ref_conditions) %>%
      column_to_rownames("sample")
    
    dds_counts <- all_counts %>%
      select(protein_id, all_of(rownames(dds_labels))) %>%
      column_to_rownames("protein_id") %>%
      as.matrix()
    
    # Factorize and relevel
    dds_labels <- dds_labels %>%
      mutate(
        strain = factor(strain),
        condition = factor(condition),
        strain_condition = factor(strain_condition)
      )
    
    dds_labels$strain <- relevel(dds_labels$strain, ref = strain_ref)
    dds_labels$condition <- relevel(dds_labels$condition, ref = cond_ref)
    dds_labels$strain_condition <- relevel(dds_labels$strain_condition, ref = strain_cond_ref)
    
    # Choose appropriate design
    design_formula <- if (length(strain_group) == 1) {
      ~ condition
    } else {
      ~ strain_condition  # captures interaction
    }
    
    # Create DESeq object
    dds <- DESeqDataSetFromMatrix(
      countData = round(dds_counts),
      colData = dds_labels,
      design = design_formula
    )
    
    # Run DESeq2
    dds <- DESeq(dds)
    
    # Save DESeq2 object
    save_path <- file.path(dds_dir, paste0("dds_", strain_set_name, "_ref_", strain_ref, "_", cond_ref, ".RData"))
    save(dds, file = save_path)
    
    
    # Save normalized + vst counts ONCE per (strain, strain_set)
    save_key <- paste0(strain_ref, "_", strain_set_name)
    if (is.null(vst_norm_saved[[save_key]])) {
      norm_counts <- counts(dds, normalized = TRUE)
      write.csv(norm_counts, file = file.path(dds_dir, paste0("norm_counts_", save_key, ".csv")))
      save(norm_counts, file = file.path(dds_dir, paste0("norm_counts_", save_key, ".RData")))
      
      vsd <- vst(dds, blind = FALSE)
      write.csv(assay(vsd), file = file.path(dds_dir, paste0("vst_", save_key, ".csv")))
      save(vsd, file = file.path(dds_dir, paste0("vst_", save_key, ".RData")))
      
      rld <- rlog(dds, blind = FALSE)
      write.csv(assay(rld), file = file.path(dds_dir, paste0("rld_", save_key, ".csv")))
      save(rld, file = file.path(dds_dir, paste0("rld_", save_key, ".RData")))
      vst_norm_saved[[save_key]] <- TRUE
    }
  }
}

message("\n All DESeq2 models and transformed counts saved.")