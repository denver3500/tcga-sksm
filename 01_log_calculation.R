###############################################################################
# Script: Log TPM Calculation and Stratification
# Description: This script processes TCGA-SKCM transcriptomic data by applying
#              log2(TPM + 1) transformation, stratifying samples based on
#              ENSG00000115419 (GLS) expression into High/Low groups, generating
#              stratification plots, and saving expression data and stratification
#              results as CSV files for Combined, Primary, and Metastatic analyses.
###############################################################################

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(SummarizedExperiment)
library(stringr)
library(readr)

# Define sample types for different analyses
sample_types_list <- list(
  "Combined" = c("Primary Tumor", "Metastatic"),
  "Primary" = c("Primary Tumor"),
  "Metastatic" = c("Metastatic")
)

# Create output directories for each analysis type
for (analysis_name in names(sample_types_list)) {
  dir.create(
    paste0("Transcriptomics-TCGA/LogTPM_plots/TCGA-SKCM-", analysis_name, "/Expression"),
    showWarnings = FALSE,
    recursive = TRUE
  )
}

# Define paths
RAW_DATA_DIR <- "raw_data"
tcga_file <- file.path(RAW_DATA_DIR, "TCGA-SKCM_transcriptomic_exp.rds")

message("Processing TCGA-SKCM dataset")

# Function to process TCGA-SKCM RDS file for a given sample type
process_tcga_rds <- function(file_path, sample_types, analysis_name) {
  message("Processing: TCGA-SKCM - ", analysis_name)

  tryCatch({
    # Load RDS file
    data <- readRDS(file_path)

    if (!inherits(data, "SummarizedExperiment")) {
      message("  Skipping - not a SummarizedExperiment object")
      return(NULL)
    }

    # Extract TPM matrix and metadata
    tpm <- assay(data, "tpm_unstrand")
    sample_info <- colData(data)
    gene_info <- rowData(data)

    # Get gene symbols
    if ("gene_name" %in% colnames(gene_info)) {
      gene_symbols <- gene_info$gene_name
      names(gene_symbols) <- rownames(tpm)
    } else {
      gene_symbols <- rownames(tpm)
      names(gene_symbols) <- rownames(tpm)
    }

    # Keep all genes (no filtering applied)
    tpm_filtered <- tpm
    message("  Found ", nrow(tpm_filtered), " genes, ", ncol(tpm_filtered), " samples")

    # Identify tumor samples based on sample types
    if ("sample_type" %in% colnames(sample_info)) {
      tumor_samples <- colnames(tpm_filtered)[sample_info$sample_type %in% sample_types]
    } else if ("shortLetterCode" %in% colnames(sample_info)) {
      code_map <- c("Primary Tumor" = "TP", "Metastatic" = "TM")
      codes <- code_map[sample_types]
      tumor_samples <- colnames(tpm_filtered)[sample_info$shortLetterCode %in% codes]
    } else {
      tumor_samples <- colnames(tpm_filtered)
    }

    if (length(tumor_samples) == 0) {
      message("  No tumor samples found for ", analysis_name)
      return(NULL)
    }

    # Filter to selected tumor samples
    tpm_filtered <- tpm_filtered[, tumor_samples, drop = FALSE]
    message("  Using ", ncol(tpm_filtered), " ", analysis_name, " samples")

    # Apply log2(TPM + 1) transformation
    log_tpm_data <- log2(tpm_filtered + 1)

    # Stratify patients based on ENSG00000115419 expression
    gene_pattern <- "^ENSG00000115419"
    gene_matches <- grep(gene_pattern, rownames(log_tpm_data))

    if (length(gene_matches) == 0) {
      message("  Gene ENSG00000115419 not found")
      return(NULL)
    }

    gene_expr <- log_tpm_data[gene_matches[1], ]
    q <- quantile(gene_expr, probs = c(0.25, 0.75))
    low_samples <- names(gene_expr)[gene_expr <= q[1]]
    high_samples <- names(gene_expr)[gene_expr >= q[2]]
    middle_samples <- names(gene_expr)[gene_expr > q[1] & gene_expr < q[2]]

    message(
      "Stratified: ", length(high_samples), " high (>=Q3), ",
      length(middle_samples), " middle (Q1-Q3), ",
      length(low_samples), " low (<=Q1) expression samples"
    )

    # Map gene IDs to symbols for plotting
    plot_gene_names <- gene_symbols[rownames(log_tpm_data)]
    plot_gene_names <- plot_gene_names[!is.na(plot_gene_names)]

    # Prepare data for stratification plot
    plot_data <- data.frame(
      sample = names(gene_expr),
      expression = gene_expr,
      group = ifelse(
        names(gene_expr) %in% high_samples, "High",
        ifelse(names(gene_expr) %in% low_samples, "Low", "Middle")
      )
    )

    # Define output directory
    output_dir <- paste0("Transcriptomics-TCGA/LogTPM_plots/TCGA-SKCM-", analysis_name, "/Expression")
    dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

    # Generate and save stratification plot
    pdf_file <- file.path(output_dir, paste0("TCGA-SKCM-", analysis_name, "_ENSG00000115419_logTPM_stratification.pdf"))

    pdf(pdf_file, width = 8, height = 6)

    p <- ggplot(plot_data, aes(x = group, y = expression, fill = group)) +
      geom_violin(alpha = 0.7) +
      geom_jitter(width = 0.2, size = 0.8, alpha = 0.5) +
      theme_bw() +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      labs(
        title = paste("TCGA-SKCM:", analysis_name, "- ENSG00000115419 Expression Stratification"),
        x = "Expression Group",
        y = "log2(TPM + 1)",
        subtitle = paste(
          "High (>=Q3):", length(high_samples), "samples, Middle (Q1-Q3):",
          length(middle_samples), "samples, Low (<=Q1):", length(low_samples), "samples"
        )
      ) +
      scale_fill_manual(values = c("High" = "red", "Low" = "blue"))

    print(p)
    dev.off()
    message("  Saved: ", pdf_file)

    # Save log TPM expression data with gene symbols
    log_tpm_with_symbols <- as.data.frame(log_tpm_data)
    log_tpm_with_symbols$gene_symbol <- plot_gene_names[rownames(log_tpm_data)]
    log_tpm_with_symbols <- log_tpm_with_symbols[!is.na(log_tpm_with_symbols$gene_symbol), ]

    write.csv(
      log_tpm_with_symbols,
      file.path(output_dir, paste0("TCGA-SKCM-", analysis_name, "_logTPM_expression.csv"))
    )

    # Save stratification information
    stratification <- data.frame(
      sample = plot_data$sample,
      ENSG00000115419_logTPM = plot_data$expression,
      group = plot_data$group
    )

    write.csv(
      stratification,
      file.path(output_dir, paste0("TCGA-SKCM-", analysis_name, "_ENSG00000115419_logTPM_stratification.csv"))
    )

    message("  Completed TCGA-SKCM - ", analysis_name)

  }, error = function(e) {
    message("  Error processing TCGA-SKCM - ", analysis_name, ": ", e$message)
  })
}

# Process the SKCM RDS file for each sample type group
for (analysis_name in names(sample_types_list)) {
  sample_types <- sample_types_list[[analysis_name]]
  process_tcga_rds(tcga_file, sample_types, analysis_name)
}

message("\nProcessing complete!")
message("Results saved in TCGA-SKCM-Combined, TCGA-SKCM-Primary, and TCGA-SKCM-Metastatic folders under LogTPM_plots")