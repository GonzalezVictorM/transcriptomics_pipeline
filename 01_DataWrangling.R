rm(list=ls())

# Load necessary library
library(tidyverse)

# Define strains of interest
strains <- c("CBS_464", "ace3_13-5")

# Function to parse GFF3 and extract gene annotations
parse_gff_genes <- function(gff_file) {
  # Read GFF, skipping comments
  gff <- read_tsv(
    gff_file,
    comment = "#",
    col_names = FALSE,
    col_types = cols(.default = "c"),
    progress = FALSE
  )
  
  # Filter only 'gene' features
  genes <- gff %>%
    filter(X3 == "gene") %>%
    select(seqid = X1, source = X2, type = X3, start = X4, end = X5, score = X6, strand = X7, phase = X8, attributes = X9)
  
  # Helper to extract values from attributes string
  extract_attr <- function(attr_str, key) {
    result <- str_match(attr_str, paste0(key, "=([^;]+)"))[,2]
    if (any(is.na(result))) message(sprintf("Warning: %d entries lack '%s' in attributes", sum(is.na(result)), key))
    result
  }
  
  # Extract specific attributes into columns
  genes <- genes %>%
    mutate(
      gene_id = extract_attr(attributes, "ID"),
      name = extract_attr(attributes, "Name"),
      product_name = extract_attr(attributes, "product_name"),
      protein_id = extract_attr(attributes, "proteinId"),
      transcript_id = extract_attr(attributes, "transcriptId")
    )
  
  # Keep only useful columns
  genes %>%
    select(gene_id, name, product_name, protein_id, transcript_id)
}

# Read the gene annotations
gene_annotations <- parse_gff_genes("00_RawData/Dicsqu464_2_GeneModels_FilteredModels2.gff3")
head(gene_annotations)

# Create output directory if it doesn't exist
output_dir <- "01_TidyData"
if (!dir.exists(output_dir)) dir.create(output_dir)

# Function to load and export data for a single strain
load_and_save_strain_data <- function(strain, gene_annotations, out_dir) {
  base_path <- file.path("00_RawData", strain)
  
  label_file <- file.path(base_path, "library_sample_match.txt")
  counts_file <- file.path(base_path, "counts.txt")
  fpkm_file <- file.path(base_path, "fpkm_counts.txt")
  
  # Checks for file existence
  files <- c(
    label = label_file,
    counts = counts_file,
    fpkm = fpkm_file
  )
  
  for (f in names(files)) {
    if (!file.exists(files[f])) stop(sprintf("File %s does not exist for strain %s", files[f], strain))
  }
  
  # Read label
  labels <- read.table(label_file, header = TRUE, sep = "\t", row.names = 1, stringsAsFactors = TRUE) %>%
    data.frame() %>%
    as_tibble(rownames = "sample")
  
  # Read and join counts
  counts <- read.table(counts_file, header = TRUE, sep = "\t") %>%
    as_tibble() %>%
    left_join(gene_annotations, by = c("Geneid" = "gene_id")) %>%
    select(-Geneid,-name,-product_name,-transcript_id) %>%
    select(protein_id, everything())
  
  # Read and join FPKMs
  fpkms <- read.table(fpkm_file, header = TRUE, sep = "\t") %>%
    as_tibble() %>%
    left_join(gene_annotations, by = c("Geneid" = "gene_id")) %>%
    select(-Geneid,-name,-product_name,-transcript_id) %>%
    select(protein_id, everything())
  
  # Save individual CSVs
  write_csv(labels, file.path(out_dir, paste0("labels_", strain, ".csv")))
  write_csv(counts, file.path(out_dir, paste0("counts_", strain, ".csv")))
  write_csv(fpkms,  file.path(out_dir, paste0("fpkms_", strain, ".csv")))
  
  list(labels = labels, counts = counts, fpkms = fpkms)
}

# Load all data and save individual files
strain_data <- set_names(strains) %>%
  map(~ load_and_save_strain_data(.x, gene_annotations, output_dir))

# Combine metadata
all_labels <- map_dfr(strain_data, "labels")

# Combine count and FPKM matrices by protein_id
all_counts <- purrr::reduce(map(strain_data, "counts"), full_join, by = "protein_id")
all_fpkms  <- purrr::reduce(map(strain_data, "fpkms"),  full_join, by = "protein_id")

# Save combined CSVs
write_csv(all_labels, file.path(output_dir, "all_strains_labels.csv"))
write_csv(all_counts, file.path(output_dir, "all_strains_counts.csv"))
write_csv(all_fpkms,  file.path(output_dir, "all_strains_fpkms.csv"))

