# ======================================================
# Spatial Transcriptomics: QC & Normalization
# Omics Core Standardized Pipeline
# ======================================================

# Usage:
# Rscript 02_qc_normalization.R <spaceranger_output_dir> <sample_id>

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 02_qc_normalization.R <spaceranger_output_dir> <sample_id>")
}

base_dir <- args[1]
sample_id <- args[2]

input_dir <- file.path(base_dir, sample_id, "outs")
output_dir <- file.path(base_dir, sample_id, "seurat_spatial")

dir.create(output_dir, showWarnings = FALSE)

cat("Loading spatial data for:", sample_id, "\n")

# ------------------------------------------------------
# Load Spatial Data
# ------------------------------------------------------

spatial_obj <- Load10X_Spatial(
  data.dir = input_dir,
  filename = "filtered_feature_bc_matrix.h5"
)

spatial_obj$orig.ident <- sample_id

# ------------------------------------------------------
# QC Metrics
# ------------------------------------------------------

spatial_obj[["percent.mt"]] <- PercentageFeatureSet(
  spatial_obj,
  pattern = "^MT-"
)

# ------------------------------------------------------
# Basic Filtering (Conservative Defaults)
# ------------------------------------------------------

spatial_obj <- subset(
  spatial_obj,
  subset = nCount_Spatial > 500 &
           percent.mt < 20
)

# ------------------------------------------------------
# Normalization & Feature Selection
# ------------------------------------------------------

spatial_obj <- SCTransform(spatial_obj, assay = "Spatial", verbose = FALSE)

# ------------------------------------------------------
# Save Object
# ------------------------------------------------------

saveRDS(
  spatial_obj,
  file = file.path(output_dir, paste0(sample_id, "_spatial_normalized.rds"))
)

# ------------------------------------------------------
# QC Plots
# ------------------------------------------------------

pdf(file.path(output_dir, paste0(sample_id, "_spatial_QC.pdf")))

SpatialFeaturePlot(spatial_obj, features = "nCount_Spatial")
SpatialFeaturePlot(spatial_obj, features = "percent.mt")

dev.off()

cat("Spatial QC & normalization complete.\n")
