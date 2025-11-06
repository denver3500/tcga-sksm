library(tidyverse)

# Read the tx2gene file
tx2gene <- read_tsv("/media/windows/BIOINFORMATICS/tcga-alexa/Cell_line/Melanoma/Nextflow/output/star_salmon/tx2gene.tsv", col_names = c("transcript_id", "gene_id", "gene_name"))

# Create unique gene_id to gene_name mapping
gene_map <- tx2gene %>%
  select(gene_id, gene_name) %>%
  distinct()

# Read the DEGs file
degs <- as.data.frame(read_csv("/media/windows/BIOINFORMATICS/tcga-alexa/Cell_line/Melanoma/Analysis/DESeq2_results.csv"))

# Set the first column as rownames
rownames(degs) <- degs[,1]
degs <- degs[, -1]

# Map Ensembl IDs to gene names for rownames
gene_names <- gene_map$gene_name[match(rownames(degs), gene_map$gene_id)]

# Handle missing mappings
gene_names[is.na(gene_names)] <- rownames(degs)[is.na(gene_names)]

# Handle duplicates: keep Ensembl ID for duplicate gene names
dup <- duplicated(gene_names) | duplicated(gene_names, fromLast = TRUE)
gene_names[dup] <- rownames(degs)[dup]

# Set new rownames
rownames(degs) <- gene_names

# Add gene column for saving
degs$gene <- rownames(degs)

# Reorder columns to put gene first
degs <- degs %>% select(gene, everything())

# Save the updated file
write_csv(degs, "/media/windows/BIOINFORMATICS/tcga-alexa/Cell_line/Melanoma/Analysis/DESeq2_results_with_gene_names.csv")

# Print summary
cat("Converted rownames from Ensembl IDs to gene names.\n")
cat("Genes without names:", sum(is.na(match(rownames(degs), gene_map$gene_id))), "\n")