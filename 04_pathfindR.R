###############################################################################
# Script: Pathway Enrichment Analysis with pathfindR
# Description: This script performs pathway enrichment analysis using pathfindR
#              on differential gene expression results from TCGA-SKCM data. It
#              analyzes GO-BP, GO-CC, and GO-MF gene sets for Combined, Primary,
#              and Metastatic analyses, generating HTML reports with enriched terms.
###############################################################################

# Load required libraries
library(pathfindR)
library(dplyr)
library(readr)

# Define gene sets for enrichment analysis
gene_sets <- c("GO-BP", "GO-CC", "GO-MF")

# Define analysis types to process
analysis_types <- c("Combined", "Primary", "Metastatic")

# Function to run pathfindR analysis for a given analysis type
run_pathfindr_analysis <- function(analysis_type) {
  message("Running pathfindR for: TCGA-SKCM - ", analysis_type)

  # Construct path to DEG CSV file
  csv_file <- file.path(
    "Transcriptomics-TCGA/LogTPM_plots",
    paste0("TCGA-SKCM-", analysis_type),
    "Expression",
    paste0("TCGA-SKCM-", analysis_type, "_DEG_High_vs_Low.csv")
  )

  if (!file.exists(csv_file)) {
    message("DEG CSV not found: ", csv_file)
    return(NULL)
  }

  # Load DEG results
  deg_data <- read.csv(csv_file)

  # Prepare input data frame for pathfindR (requires Gene.symbol, logFC, adj.P.Val)
  input_df <- deg_data %>%
    select(Gene.symbol = gene_symbol, logFC = log2FoldChange, adj.P.Val = padj) %>%
    filter(!is.na(Gene.symbol) & !is.na(logFC) & !is.na(adj.P.Val))

  if (nrow(input_df) == 0) {
    message("No valid data for pathfindR in ", analysis_type)
    return(NULL)
  }

  # Define output directory
  output_dir <- file.path(
    "Transcriptomics-TCGA/LogTPM_plots",
    paste0("TCGA-SKCM-", analysis_type),
    "Expression"
  )

  # Run pathfindR for each gene set
  for (gene_set in gene_sets) {
    message("  Processing gene set: ", gene_set)
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

  message("Completed pathfindR for ", analysis_type)
}

# Run pathfindR analysis for each analysis type
for (type in analysis_types) {
  run_pathfindr_analysis(type)
}

message("pathfindR analysis complete!")