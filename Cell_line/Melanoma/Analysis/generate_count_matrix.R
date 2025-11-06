library(tidyverse)
library(DESeq2)
library(tximport)

# Set working directory
setwd("/media/windows/BIOINFORMATICS/tcga-alexa/Cell_line/Melanoma")

# Load metadata
meta <- read_csv("sample_sheet.csv")

# Directory where salmon gene quant files are
dir <- "output/star_salmon"

# List all directories containing quant.genes.sf files using the sample column of metadata
files <- file.path(dir, meta$sample, "quant.genes.sf")

# Name the file list with the samplenames
names(files) <- meta$sample

# Check if all files exist
if (!all(file.exists(files))) {
  stop("Some quant.genes.sf files do not exist")
}

# Since tximport may not work directly with gene-level, let's do it manually
# Read all quant files and combine NumReads

# Read the first file to get gene names
first_file <- files[1]
data <- read_tsv(first_file, col_types = cols()) %>%
  select(Name, NumReads) %>%
  setNames(c("Name", names(files)[1]))

# Read the rest and join
for (i in 2:length(files)) {
  temp <- read_tsv(files[i], col_types = cols()) %>%
    select(Name, NumReads) %>%
    setNames(c("Name", names(files)[i]))
  data <- full_join(data, temp, by = "Name")
}

# Set row names
data <- data %>% column_to_rownames("Name")

# Round to nearest whole number
data <- round(data)

# Remove rows with all zeros
keep <- rowSums(data, na.rm = TRUE) > 0
data <- data[keep, ]

# Save the count matrix
write.csv(data, "count_matrix.csv", row.names = TRUE)