###############################################################################
# Script: Download TCGA Transcriptomic Data
# Description: This script downloads transcriptomic expression data from TCGA
#              for SKCM using TCGAbiolinks. It queries, downloads,
#              and prepares the data and then saving it as RDS file.
###############################################################################

# Load required libraries
library(SummarizedExperiment)
library(TCGAbiolinks)

# Define projects to process (only SKCM here)
project_ids <- c("TCGA-SKCM")

# Loop through each project to query, download, and prepare data
for (project_id in project_ids) {
  message(paste("Processing project:", project_id))

  # Query TCGA data
  query <- tryCatch({
    GDCquery(
      project = project_id,
      data.category = "Transcriptome Profiling",
      data.type = "Gene Expression Quantification",
      workflow.type = "STAR - Counts",
      sample.type = c("Primary Tumor", "Solid Tissue Normal", "Metastatic")
    )
  }, error = function(e) {
    message(paste("Error creating query for project:", project_id))
    message(e)
    return(NULL)
  })

  if (!is.null(query)) {
    # Download the queried data
    tryCatch({
      GDCdownload(
        directory = "raw_data/GDCdata",
        query = query,
        files.per.chunk = 100
      )
    }, error = function(e) {
      message(paste("Error downloading data for project:", project_id))
      message(e)
    })

    # Prepare and save the data
    tryCatch({
      transcriptomic_exp <- GDCprepare(
        directory = "raw_data/GDCdata",
        query = query
      )

      # Save the prepared data as summarized experiment in RDS file
      saveRDS(
        transcriptomic_exp,
        file = file.path("raw_data", paste0(project_id, "_transcriptomic_exp.rds"))
      )
    }, error = function(e) {
      message(paste("Error preparing data for project:", project_id))
      message(e)
    })
  }
}