rm(list = ls())

# Load required libraries
library(tidyverse)
library(DESeq2)
library(apeglm)

# Set analysis parameters
padj.cutoff <- 0.05
fc.cutoff <- 2.5

# Define paths
dds_dir <- "03_DESeqObjects"
deg_dir <- "04_DEGanalysis"
if (!dir.exists(deg_dir)) dir.create(deg_dir)

# DEG classification function
DEG <- function(x) {
  mode(x) <- "numeric"
  
  if (is.na(x[1]) | is.na(x[2])) {
    return("Lowly expressed")
  } else if (x[1] >= fc.cutoff & x[2] < padj.cutoff) {
    return("Overexpressed")
  } else if (x[1] <= (1 / fc.cutoff) & x[2] < padj.cutoff) {
    return("Underexpressed")
  } else {
    return("Non-significant")
  }
}

# List all DESeq2 model files
dds_files <- list.files(dds_dir, pattern = "^dds_.*\\.RData$", full.names = TRUE)

#dds_file <- dds_files[20] # used for debugging
# Loop through each model file
for (dds_file in dds_files) {
  message("🔍 Processing: ", basename(dds_file))
  
  # Load the DESeq2 object
  load(dds_file)
  
  if (!exists("dds") || class(dds)[1] != "DESeqDataSet") {
    warning("Skipping file — no valid 'dds' object found in ", dds_file)
    next
  }
  
  # Extract model coefficients (excluding intercept)
  coef_names <- resultsNames(dds)[-1]
  
  # Build simplified contrast labels
  contrast_labels <- sapply(coef_names, function(x) {
    gsub("strain_condition_", "", x)
  }, USE.NAMES = TRUE)
  
  # Prepare result table
  resLFC_table <- data.frame()
  
  #coef <- coef_names[16]# used for debugging
  
  for (coef in coef_names) {
    # Get shrunken results
    res <- results(dds, name = coef, alpha = padj.cutoff)
    resLFC <- lfcShrink(dds, coef = coef, res = res, type = "apeglm") %>%
      as.data.frame() %>%
      mutate(
        FoldChange = 2^log2FoldChange,
        cSource = contrast_labels[coef]
      ) %>%
      rownames_to_column("ProteinID")
    
    resLFC_table <- bind_rows(resLFC_table, resLFC)
  }
  
  # Classify genes based on fold change and adjusted p-value
  resLFC_table$Expression <- apply(
    resLFC_table[, c("FoldChange", "padj")],
    1,
    DEG
  )
  
  # Output file name
  out_file <- file.path(
    deg_dir,
    paste0(tools::file_path_sans_ext(basename(dds_file)), "_DEG_p", padj.cutoff, "_fc", fc.cutoff, ".csv")
  )
  
  # Save to CSV
  write_csv(resLFC_table, out_file)
}

message("✅ DEG analysis complete. Results saved to ", deg_dir)
