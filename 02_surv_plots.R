###############################################################################
# Script: Survival Analysis Plots
# Description: This script generates Kaplan-Meier survival plots for TCGA-SKCM
#              patients stratified by ENSG00000115419 (GLS) expression levels (High/Low)
#              and combined with AJCC pathologic stages (Early/Late). It produces
#              separate plots for expression-only and combined stratification for
#              Combined, Primary, and Metastatic analyses.
###############################################################################

# Load required libraries
library(TCGAbiolinks)
library(dplyr)

# Retrieve clinical data for TCGA-SKCM
clin.skcm <- GDCquery_clinic("TCGA-SKCM", "clinical")

# Define analysis types to process
analysis_types <- c("Combined", "Primary", "Metastatic")

# Function to generate survival plots for a given analysis type
plot_survival <- function(analysis_type) {
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

  # Filter to High and Low expression groups only (exclude Middle)
  strat_data <- strat_data[strat_data$group %in% c("High", "Low"), ]

  # Extract patient ID from sample ID (first 12 characters)
  strat_data$patient <- substr(strat_data$sample, 1, 12)

  # Merge stratification data with clinical data
  clin_merged <- merge(
    clin.skcm,
    strat_data[, c("patient", "group")],
    by.x = "submitter_id",
    by.y = "patient",
    all.x = FALSE
  )

  if (nrow(clin_merged) == 0) {
    message("No matching clinical data for ", analysis_type)
    return(NULL)
  }

  # Generate survival plot for expression groups only
  pdf_file <- file.path(
    "Transcriptomics-TCGA/LogTPM_plots",
    paste0("TCGA-SKCM-", analysis_type),
    "Expression",
    paste0("TCGA-SKCM-", analysis_type, "_survival_plot.pdf")
  )

  TCGAanalyze_survival(
    data = clin_merged,
    clusterCol = "group",
    main = paste("TCGA SKCM\n", analysis_type, "- ENSG00000115419 Stratification"),
    height = 10,
    width = 10,
    filename = pdf_file
  )

  message("Saved survival plot (expression): ", pdf_file)

  # Define AJCC stage groupings
  early_stages <- c("Stage 0", "Stage I", "Stage IA", "Stage IB", "Stage II", "Stage IIA", "Stage IIB")
  late_stages <- c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "Stage IV")

  # Categorize patients into Early/Late stages
  clin_merged$stage_group <- ifelse(
    clin_merged$ajcc_pathologic_stage %in% early_stages, "Early",
    ifelse(clin_merged$ajcc_pathologic_stage %in% late_stages, "Late", NA)
  )

  # Remove patients with missing stage information
  clin_merged <- clin_merged[!is.na(clin_merged$stage_group), ]

  if (nrow(clin_merged) == 0) {
    message("No valid stage data for ", analysis_type)
    return(NULL)
  }

  # Create combined stratification groups (e.g., High_Early, Low_Late)
  clin_merged$combined_group <- paste(clin_merged$group, clin_merged$stage_group, sep = "_")

  # Filter to combined groups with at least 5 patients
  group_counts <- table(clin_merged$combined_group)
  valid_groups <- names(group_counts)[group_counts >= 5]
  clin_merged <- clin_merged[clin_merged$combined_group %in% valid_groups, ]

  if (nrow(clin_merged) == 0) {
    message("No valid combined groups with sufficient data for ", analysis_type)
    return(NULL)
  }

  # Generate survival plot for combined stratification
  pdf_file_combined <- file.path(
    "Transcriptomics-TCGA/LogTPM_plots",
    paste0("TCGA-SKCM-", analysis_type),
    "Expression",
    paste0("TCGA-SKCM-", analysis_type, "_survival_plot_combined.pdf")
  )

  TCGAanalyze_survival(
    data = clin_merged,
    clusterCol = "combined_group",
    main = paste("TCGA SKCM\n", analysis_type, "- Combined Stratification"),
    height = 10,
    width = 10,
    conf.int = FALSE,
    filename = pdf_file_combined
  )

  message("Saved survival plot (combined): ", pdf_file_combined)
}

# Generate survival plots for each analysis type
for (type in analysis_types) {
  plot_survival(type)
}