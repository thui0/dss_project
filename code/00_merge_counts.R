# source this script in "DSS_Combined_Project.Rproj" directory
# for more details see accompanying .Rmd file "merge_counts_tables_nb.Rmd"

# libraries needed for this script
library(dplyr)

# read counts tables
x <- read.csv(file = file.path("data", "00_batch_counts_tables", "batch_01_02_merged_counts_table.csv"))
y <- read.csv(file = file.path("data", "00_batch_counts_tables", "batch_03_counts_table.csv"))

# arrange by "gene_id"
x <- x %>% arrange(gene_id)
y <- y %>% arrange(gene_id)

# merge the counts tables
z <- merge(x = x,
           y = y,
           by.x = "gene_id",
           by.y = "gene_id",
           all = TRUE,
)

# remove rows with NAs
z <- z %>% na.omit()

# rearrange columns
z.new <- z %>% 
  relocate(MP2, .after = Z3) %>%
  relocate(MP1, .after = MS1)

# save merged counts table

## save as RDS file for DGE analysis input
# saveRDS(z.new, file = file.path("data", "merged_counts_table.rds"))

## save as CSV file for DGE analysis input
write.csv(z.new, file = file.path("data", "merged_counts_table.csv"), row.names = FALSE)

## save as tab-delimited TXT file in sub-directory for DGE analysis input
# write.table(z.new, file = file.path("data", "merged_counts_table.txt"), row.names = FALSE)