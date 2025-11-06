library(tidyverse)

# Read the tx2gene file
tx2gene <- read_tsv("output/star_salmon/tx2gene.tsv", col_names = c("transcript_id", "gene_id", "gene_name"))

# Create unique gene_id to gene_name mapping
gene_map <- tx2gene %>%
  select(gene_id, gene_name) %>%
  distinct()

# Read the DEGs file
degs <- read_csv("all_DEGs.csv")

# Merge to add gene names
degs_with_names <- degs %>%
  left_join(gene_map, by = c("gene" = "gene_id"))

# Save the updated file
write_csv(degs_with_names, "all_DEGs_with_names.csv")

# Print summary
cat("Added gene names to", nrow(degs_with_names), "genes.\n")
cat("Genes without names:", sum(is.na(degs_with_names$gene_name)), "\n")