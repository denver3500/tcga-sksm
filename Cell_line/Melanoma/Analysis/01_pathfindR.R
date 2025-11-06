###############################################################################
# Script: Pathway Enrichment Analysis with pathfindR
# Description: This script performs pathway enrichment analysis using pathfindR
#              on differential gene expression results from cell line data.
###############################################################################

# Load required libraries
library(pathfindR)
library(dplyr)
library(readr)

# Define gene sets for enrichment analysis
gene_sets <- c("GO-BP", "GO-CC", "GO-MF")

# Path to DEG CSV file
csv_file <- "/media/windows/BIOINFORMATICS/tcga-alexa/Cell_line/Melanoma/Analysis/DESeq2_results_with_gene_names.csv"

if (!file.exists(csv_file)) {
  stop("DEG CSV not found: ", csv_file)
}

# Load DEG results
deg_data <- read_csv(csv_file)

# Prepare input data frame for pathfindR (requires Gene.symbol, logFC, adj.P.Val)
input_df <- deg_data %>%
  select(Gene.symbol = gene, logFC = log2FoldChange, adj.P.Val = padj) %>%
  filter(complete.cases(.))

# Ensure correct data types
input_df$Gene.symbol <- as.character(input_df$Gene.symbol)
input_df$logFC <- as.numeric(input_df$logFC)
input_df$adj.P.Val <- as.numeric(input_df$adj.P.Val)

input_df <- as.data.frame(input_df)

# Debug
print("Class of adj.P.Val:")
print(class(input_df$adj.P.Val))
print("Any NA in adj.P.Val:")
print(any(is.na(input_df$adj.P.Val)))
print("Head of input_df:")
print(head(input_df))

if (nrow(input_df) == 0) {
  stop("No valid data for pathfindR")
}

# Define output directory
output_dir <- "pathfindR_output"

# Run pathfindR for each gene set
for (gene_set in gene_sets) {
  message("Processing gene set: ", gene_set)
  run_pathfindR(
    input = input_df,
    gene_sets = gene_set,
    pin_name = "Biogrid",
    min_gset_size = 10,
    max_gset_size = 300,
    p_val_threshold = 0.05,
    output_dir = output_dir
  )
}

message("pathfindR analysis complete!")