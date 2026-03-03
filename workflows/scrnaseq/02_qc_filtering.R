# ======================================================
# Single-Cell RNA-seq: QC Filtering (Seurat)
# Omics Core Standardized Pipeline
# ======================================================

# Usage:
# Rscript 02_qc_filtering.R <cellranger_output_dir> <sample_id>

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 02_qc_filtering.R <cellranger_output_dir> <sample_id>")
}

base_dir <- args[1]
sample_id <- args[2]

input_dir <- file.path(base_dir, sample_id, "outs", "filtered_feature_bc_matrix")
output_dir <- file.path(base_dir, sample_id, "seurat_objects")

dir.create(output_dir, showWarnings = FALSE)

# ------------------------------------------------------
# Load Data
# ------------------------------------------------------

cat("Loading data for sample:", sample_id, "\n")

counts <- Read10X(data.dir = input_dir)
seurat_obj <- CreateSeuratObject(counts = counts, project = sample_id)

# ------------------------------------------------------
# Calculate QC Metrics
# ------------------------------------------------------

seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

# ------------------------------------------------------
# Basic Filtering (Core defaults)
# ------------------------------------------------------

seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > 200 &
           nFeature_RNA < 6000 &
           percent.mt < 15
)

# ------------------------------------------------------
# Save Filtered Object
# ------------------------------------------------------

saveRDS(seurat_obj,
        file = file.path(output_dir, paste0(sample_id, "_qc_filtered.rds")))

# ------------------------------------------------------
# QC Plots
# ------------------------------------------------------

pdf(file.path(output_dir, paste0(sample_id, "_QC_violin.pdf")))
VlnPlot(seurat_obj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        ncol = 3)
dev.off()

cat("QC filtering complete for sample:", sample_id, "\n")
