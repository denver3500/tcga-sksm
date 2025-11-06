library(DESeq2)
library(tidyverse)
library(pheatmap)
library(sva)

# Read count matrix
counts <- read.csv("/media/windows/BIOINFORMATICS/tcga-alexa/Cell_line/Melanoma/count_matrix.csv", row.names = 1)

# Read sample metadata
meta <- read_csv("/media/windows/BIOINFORMATICS/tcga-alexa/Cell_line/Melanoma/sample_sheet.csv")

# Create condition column based on sample names
meta <- meta %>%
  mutate(condition = ifelse(grepl("SCR", sample), "SCR", "GLS_KO"))

# Create batch factor
batch <- factor(ifelse(grepl("exp1|exp2", meta$sample), 1, 2))

# Adjust counts for batch effects using ComBat_seq
adjusted_counts <- sva::ComBat_seq(counts, batch = batch, group = NULL)

# Ensure samples are in the same order
adjusted_counts <- adjusted_counts[, meta$sample]

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = adjusted_counts,
                              colData = meta,
                              design = ~ condition)

# Run DESeq
dds <- DESeq(dds)

# Get differential expression results
res <- results(dds, contrast = c("condition", "GLS_KO", "SCR"))

# Save results
write.csv(as.data.frame(res), "DESeq2_results.csv")

# Optionally, get significant genes
sig_res <- res[!is.na(res$padj) & res$padj < 0.05 & abs(res$log2FoldChange) > 1, ]
write.csv(as.data.frame(sig_res), "significant_DEGs.csv")

# Get variance stabilized counts
vst <- vst(dds, blind = TRUE)

# Compute Euclidean distance
sampleDists <- dist(t(assay(vst)), method = "euclidean")

# Convert to matrix
sampleDistMatrix <- as.matrix(sampleDists)

# Plot heatmap
ann_colors <- list(condition = c(SCR = "blue", GLS_KO = "red"))
ann_df <- meta %>% select(condition) %>% as.data.frame()
rownames(ann_df) <- meta$sample
pheatmap(sampleDistMatrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = ann_df,
         annotation_row = ann_df,
         annotation_colors = ann_colors,
         main = "Euclidean Distance Heatmap of Samples",
         filename = "euclidean_distance_heatmap.png")