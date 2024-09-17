# filter DESeq2 results
# keep only genes that have padj < 0.05

library(tidyverse)

res_path <- file.path("results", "dge_analysis", "diffexp")
out_dir <- file.path("results", "filter_dge_results")
if (!dir.exists(out_dir)) {
  dir.create(out_dir, showWarnings = FALSE)
}

## create helper function
filter_and_save <- function(rds_path, comparison_name) {
  res <- readRDS(rds_path)
  res.filtered <- data.frame(res) |> filter(padj < 0.05) |> arrange(padj)
  print(dim(res.filtered))
  write.csv(res.filtered, file = file.path(out_dir,paste0(comparison_name,"_filtered.csv")))
  sink(file = file.path(out_dir,paste0(comparison_name,"_DEG_summary.txt")))
  print(summary(res))
  sink()
}

# KODSS vs KOnoDSS
filter_and_save(rds_path = file.path(res_path, "KODSS_vs_KOnoDSS.rds"),
                comparison_name = "KODSS_vs_KOnoDSS")

# KODSS vs WTDSS
filter_and_save(rds_path = file.path(res_path, "KODSS_vs_WTDSS.rds"),
                comparison_name = "KODSS_vs_WTDSS")

# KODSS vs WTnoDSS
filter_and_save(rds_path = file.path(res_path, "KODSS_vs_WTnoDSS.rds"),
                comparison_name = "KODSS_vs_WTnoDSS")

# KOnoDSS vs WTnoDSS
filter_and_save(rds_path = file.path(res_path, "KOnoDSS_vs_WTnoDSS.rds"),
                comparison_name = "KOnoDSS_vs_WTnoDSS")

# WTDSS vs WTnoDSS
filter_and_save(rds_path = file.path(res_path, "WTDSS_vs_WTnoDSS.rds"),
                comparison_name = "WTDSS_vs_WTnoDSS")

# interaction term
filter_and_save(rds_path = file.path(res_path, "genotypeKO.treatmentDSS_interaction_term.rds"),
                comparison_name = "genotypeKO.treatmentDSS_interaction_term")

# filter norm counts
norm.counts = readRDS(file.path("results", "dge_analysis", "norm.counts.rds"))

# remove rows with row sum of 0
norm.counts.filtered = subset(norm.counts, rowSums(norm.counts[, 1:(ncol(norm.counts) - 1)]) > 0)

# remove the "mapIds" column from the dataframe
norm.counts.filtered <- norm.counts.filtered[, !names(norm.counts.filtered) %in% "mapIds"]

norm.counts.filtered <- norm.counts.filtered %>% 
  rownames_to_column(var = "gene_id") %>% 
  select(gene_id, everything())

saveRDS(norm.counts.filtered, file = file.path(out_dir, "norm.counts.filtered.rds"))
write.table(norm.counts.filtered, file = file.path(out_dir, "norm.counts.filtered.txt"), sep = "\t", quote = F, row.names = F)
