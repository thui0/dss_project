# create phenotype label CLS file for Standard GSEA

library(tidyverse)

out.dir = file.path("results", "gsea")
dir.create(out.dir, showWarnings = FALSE)

# read in metadata table
metadata = read.csv(file.path("data", "sample_metadata.csv"), header = T, row.names = 1)

# read in normalized counts table
norm_counts = read.table(file.path("results", "filter_dge_results", "norm.counts.filtered.txt"), header = T)

# Extract the sample names from the data frame column names
sample_names <- colnames(norm_counts)[-1]  # Excluding the "gene_id" column

# Filter the metadata for the samples in 'norm.counts.filtered'
filtered_metadata <- metadata %>% 
  rownames_to_column("sample") %>% 
  filter(sample %in% sample_names)

# Create the combined phenotype label
filtered_metadata$phenotype <- paste0(filtered_metadata$genotype, filtered_metadata$treatment)

# Ensure that the order of samples in the metadata matches the order in 'norm.counts'
filtered_metadata <- filtered_metadata %>%
  mutate(sample = factor(sample, levels = sample_names)) %>%
  arrange(sample)

# Check if the order actually matches
all(filtered_metadata$sample == sample_names)

# Create the .cls file content
num_samples <- nrow(filtered_metadata)
num_classes <- length(unique(filtered_metadata$phenotype))
class_names <- paste(unique(filtered_metadata$phenotype), collapse = "\t")

# Sample labels
sample_labels <- paste(filtered_metadata$phenotype, collapse = "\t")

# Write the .cls file
writeLines(c(
  paste(num_samples, num_classes, 1),
  paste0("#","\t",class_names),
  sample_labels
), file.path(out.dir, "DSS_combined_phenotype.cls"))

################################################################################

# create expression dataset GCT file for Standard GSEA


# Add a "Description" column filled with "na"
norm_counts$Description <- "na"

# Ensure that "Description" is the second column
norm_counts <- norm_counts %>% select(gene_id, Description, everything())

# Rename the "gene_id" column to "NAME"
names(norm_counts)[names(norm_counts) == "gene_id"] <- "NAME"

# Add a "Description" column filled with "na"
norm_counts$Description <- "na"

# Ensure that "Description" is the second column
norm_counts <- norm_counts %>% select(NAME, Description, everything())

# Function to create a .gct file
write_gct <- function(df, filename) {
  num_genes <- nrow(df)
  num_samples <- ncol(df) - 2  # Subtract 2 for the "NAME" and "Description" columns
  
  # Create the headers for the .gct file
  headers <- c("#1.2", paste(num_genes, num_samples))
  
  # Create the descriptive headers for the gene expression data
  desc_headers <- colnames(df)
  
  # Write the .gct file
  writeLines(headers, filename)
  write.table(rbind(desc_headers, df), file = filename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
}

# Use the function to create the .gct file
write_gct(norm_counts, file.path(out.dir, "DSS_combined_norm_counts_filtered.gct"))
