# ======================================================
# Spatial Transcriptomics: Clustering & Spatial Mapping
# Omics Core Standardized Pipeline
# ======================================================

# Usage:
# Rscript 03_clustering_spatial_mapping.R <spaceranger_output_dir> <sample_id>

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 03_clustering_spatial_mapping.R <spaceranger_output_dir> <sample_id>")
}

base_dir <- args[1]
sample_id <- args[2]

input_path <- file.path(
  base_dir,
  sample_id,
  "seurat_spatial",
  paste0(sample_id, "_spatial_normalized.rds")
)

if (!file.exists(input_path)) {
  stop("Normalized spatial object not found.")
}

output_dir <- file.path(base_dir, sample_id, "spatial_analysis")
dir.create(output_dir, showWarnings = FALSE)

cat("Loading normalized spatial object...\n")

spatial_obj <- readRDS(input_path)

# ------------------------------------------------------
# Dimensionality Reduction
# ------------------------------------------------------

spatial_obj <- RunPCA(spatial_obj, verbose = FALSE)
spatial_obj <- FindNeighbors(spatial_obj, dims = 1:20)
spatial_obj <- FindClusters(spatial_obj, resolution = 0.5)
spatial_obj <- RunUMAP(spatial_obj, dims = 1:20)

# ------------------------------------------------------
# Save Clustered Object
# ------------------------------------------------------

saveRDS(
  spatial_obj,
  file = file.path(output_dir, paste0(sample_id, "_spatial_clustered.rds"))
)

# ------------------------------------------------------
# Visualization
# ------------------------------------------------------

pdf(file.path(output_dir, paste0(sample_id, "_UMAP_clusters.pdf")))
DimPlot(spatial_obj, reduction = "umap", label = TRUE)
dev.off()

pdf(file.path(output_dir, paste0(sample_id, "_Spatial_clusters.pdf")))
SpatialDimPlot(spatial_obj, label = TRUE, label.size = 3)
dev.off()

cat("Spatial clustering & mapping complete.\n")
