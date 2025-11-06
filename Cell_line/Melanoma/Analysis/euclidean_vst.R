library(tidyverse)
library(pheatmap)

# Read the VST normalized counts
lines <- readLines("/media/windows/BIOINFORMATICS/tcga-alexa/Cell_line/Melanoma/double_check/vst_normalized_counts.txt")
header <- strsplit(lines[1], " ")[[1]]
header <- gsub("\"", "", header)
data <- read.table(textConnection(lines[-1]), sep = " ", header = FALSE)
colnames(data) <- header
rownames(data) <- as.character(data[,1])
vst_counts <- data[,-1]

# Ensure all columns are numeric
vst_counts <- as.data.frame(lapply(vst_counts, as.numeric))

# Transpose so samples are rows
vst_t <- t(vst_counts)

# Compute Euclidean distance
sampleDists <- dist(vst_t, method = "euclidean")

# Convert to matrix
sampleDistMatrix <- as.matrix(sampleDists)

# Create annotation for conditions
sample_conditions <- sub("_\\d+$", "", rownames(vst_t))
ann_df <- data.frame(condition = factor(sample_conditions, levels = c("CTRL", "TRAP1_KO", "GLS1_KO", "DKO")))
rownames(ann_df) <- rownames(vst_t)

ann_colors <- list(condition = c(CTRL = "green", TRAP1_KO = "blue", GLS1_KO = "red", DKO = "purple"))

# Plot heatmap
pheatmap(sampleDistMatrix,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         annotation_col = ann_df,
         annotation_colors = ann_colors,
         main = "Euclidean Distance Heatmap of Samples (VST Counts)",
         filename = "euclidean_distance_heatmap_vst.png")