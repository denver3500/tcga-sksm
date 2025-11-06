library(DESeq2)
library(tidyverse)

# Read count matrix
counts <- read.csv("count_matrix.csv", row.names = 1)

# Read sample metadata
meta <- read_csv("sample_sheet.csv")

# Create condition column based on sample names
meta <- meta %>%
  mutate(condition = ifelse(grepl("SCR", sample), "SCR", "GLS_KO"))

# Ensure samples are in the same order
counts <- counts[, meta$sample]

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = meta,
                              design = ~ condition)

# Run DESeq
dds <- DESeq(dds)

# Get results
res <- results(dds, contrast = c("condition", "GLS_KO", "SCR"))

# Summary
summary(res)

# Convert to dataframe and add gene column
res_df <- res %>%
  as.data.frame() %>%
  rownames_to_column("gene")

# Save all DEGs
write.csv(res_df, "all_DEGs.csv", row.names = FALSE)

# Print number of genes
cat("Total number of genes analyzed:", nrow(res_df), "\n")