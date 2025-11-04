###############################################################################
# Script: Differential Gene Expression Analysis
# Description: This script performs differential gene expression (DEG) analysis
#              using DESeq2 on TCGA-SKCM transcriptomic data. It compares High vs
#              Low expression groups of ENSG00000115419 (GLS) for Combined, Primary, and
#              Metastatic analyses, saving results as CSV files with gene symbols.
###############################################################################

# Load required libraries
library(DESeq2)
library(SummarizedExperiment)
library(dplyr)
library(readr)

# Define analysis types to process
analysis_types <- c("Combined", "Primary", "Metastatic")

# Path to the TCGA-SKCM RDS file
tcga_file <- "raw_data/TCGA-SKCM_transcriptomic_exp.rds"

# Load the transcriptomic data
data <- readRDS(tcga_file)

# Function to perform DEG analysis for a given analysis type
perform_deg <- function(analysis_type) {
  message("Performing DEG for: TCGA-SKCM - ", analysis_type)

  # Construct path to stratification CSV file
  csv_file <- file.path(
    "Transcriptomics-TCGA/LogTPM_plots",
    paste0("TCGA-SKCM-", analysis_type),
    "Expression",
    paste0("TCGA-SKCM-", analysis_type, "_ENSG00000115419_logTPM_stratification.csv")
  )

  if (!file.exists(csv_file)) {
    message("CSV file not found: ", csv_file)
    return(NULL)
  }

  # Load stratification data
  strat_data <- read.csv(csv_file)

  # Filter to High and Low expression groups only
  strat_data <- strat_data[strat_data$group %in% c("High", "Low"), ]

  if (nrow(strat_data) == 0) {
    message("No High/Low samples for ", analysis_type)
    return(NULL)
  }

  # Extract sample IDs for each group
  high_samples <- strat_data$sample[strat_data$group == "High"]
  low_samples <- strat_data$sample[strat_data$group == "Low"]

  message("  High samples: ", length(high_samples), ", Low samples: ", length(low_samples))

  # Subset data to selected samples
  selected_samples <- c(high_samples, low_samples)
  data_subset <- data[, colnames(data) %in% selected_samples]

  # Extract raw count matrix
  counts <- assay(data_subset)

  # Create column data for DESeq2
  colData <- data.frame(
    sample = colnames(counts),
    group = ifelse(colnames(counts) %in% high_samples, "High", "Low")
  )
  rownames(colData) <- colData$sample

  # Create DESeq2 dataset
  dds <- DESeqDataSetFromMatrix(
    countData = counts,
    colData = colData,
    design = ~ group
  )

  # Filter genes with low counts (at least 10 counts in at least 10% of samples)
  keep <- rowSums(counts(dds) >= 10) >= ceiling(0.1 * ncol(counts))
  dds <- dds[keep, ]

  # Run DESeq2 analysis
  dds <- DESeq(dds)

  # Extract results for High vs Low comparison
  res <- results(dds, contrast = c("group", "High", "Low"))

  # Add gene symbols to results
  gene_info <- rowData(data_subset)
  if ("gene_name" %in% colnames(gene_info)) {
    res$gene_symbol <- gene_info$gene_name[match(rownames(res), rownames(gene_info))]
  }

  # Convert results to data frame and sort by adjusted p-value
  res_df <- as.data.frame(res)
  res_df <- res_df[order(res_df$padj), ]

  # Define output directory
  output_dir <- file.path(
    "Transcriptomics-TCGA/LogTPM_plots",
    paste0("TCGA-SKCM-", analysis_type),
    "Expression"
  )

  # Save DEG results to CSV
  write.csv(
    res_df,
    file.path(output_dir, paste0("TCGA-SKCM-", analysis_type, "_DEG_High_vs_Low.csv"))
  )

  message("  Saved DEG results: ", nrow(res_df), " genes")
}

# Perform DEG analysis for each analysis type
for (type in analysis_types) {
  perform_deg(type)
}

message("DEG analysis complete!")