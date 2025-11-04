###############################################################################
# Script: Generate Enrichment Charts from pathfindR Results
# Description: This script loads enriched terms CSV files from pathfindR analysis
#              and generates enrichment charts for GO-BP, GO-CC, and GO-MF gene sets
#              across Combined, Primary, and Metastatic analyses. Charts are saved
#              as PDF files showing top enriched terms.
###############################################################################

# Load required libraries
library(pathfindR)
library(readr)
library(ggplot2)

# Define analysis types and gene sets to process
analysis_types <- c("Combined", "Primary", "Metastatic")
gene_sets <- c("GO-BP", "GO-CC", "GO-MF")

# Function to process a single enriched terms CSV and generate enrichment chart
process_csv <- function(analysis_type, gene_set) {
  message("Processing: ", analysis_type, " - ", gene_set)

  # Construct path to enriched terms CSV file
  csv_file <- paste0(
    "Transcriptomics-TCGA/LogTPM_plots/TCGA-SKCM-", analysis_type,
    "/Expression/TCGA-SKCM-", analysis_type, "_", gene_set, "_enriched_terms.csv"
  )

  if (!file.exists(csv_file)) {
    message("CSV not found: ", csv_file)
    return(NULL)
  }

  # Load enriched terms data
  output_df <- as.data.frame(read_csv(csv_file, show_col_types = FALSE))

  # Generate enrichment chart (top 15 terms)
  clustered_df <- enrichment_chart(output_df, top_terms = 15)

  # Define output PDF file path
  pdf_file <- paste0(
    "Transcriptomics-TCGA/LogTPM_plots/TCGA-SKCM-", analysis_type,
    "/Expression/TCGA-SKCM-", analysis_type, "_", gene_set, "_enrichment_chart.pdf"
  )

  # Save chart as PDF
  ggsave(pdf_file, clustered_df)

  message("Saved: ", pdf_file)
}

# Process all combinations of analysis types and gene sets
for (type in analysis_types) {
  for (gs in gene_sets) {
    process_csv(type, gs)
  }
}

message("All enrichment charts generated!")