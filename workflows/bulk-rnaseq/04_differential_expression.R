# ======================================================
# Bulk RNA-seq: Differential Expression (DESeq2)
# Omics Core Standardized Pipeline
# ======================================================

# Usage:
# Rscript 04_differential_expression.R counts.txt metadata.csv

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 04_differential_expression.R gene_counts.txt metadata.csv")
}

count_file <- args[1]
metadata_file <- args[2]

# ------------------------------------------------------
# 1. Import Data
# ------------------------------------------------------

counts <- read.delim(count_file, comment.char = "#")

# Remove annotation columns from featureCounts
counts_clean <- counts[, -c(1:6)]
rownames(counts_clean) <- counts$Geneid

metadata <- read.csv(metadata_file, row.names = 1)

# Ensure column order matches metadata
counts_clean <- counts_clean[, rownames(metadata)]

# ------------------------------------------------------
# 2. Create DESeq2 Object
# ------------------------------------------------------

dds <- DESeqDataSetFromMatrix(
  countData = counts_clean,
  colData = metadata,
  design = ~ condition
)

# Filter low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# ------------------------------------------------------
# 3. Run Differential Expression
# ------------------------------------------------------

dds <- DESeq(dds)

res <- results(dds)
res_ordered <- res[order(res$padj), ]

# ------------------------------------------------------
# 4. Export Results
# ------------------------------------------------------

write.csv(as.data.frame(res_ordered),
          file = "differential_expression_results.csv")

normalized_counts <- counts(dds, normalized = TRUE)

write.csv(as.data.frame(normalized_counts),
          file = "normalized_counts.csv")

# ------------------------------------------------------
# 5. QC: PCA Plot
# ------------------------------------------------------

vsd <- vst(dds, blind = FALSE)

pdf("PCA_plot.pdf")
plotPCA(vsd, intgroup = "condition")
dev.off()

message("Differential expression analysis complete.")
