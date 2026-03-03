# ======================================================
# Spatial Transcriptomics: Marker & Region Analysis
# Omics Core Standardized Pipeline
# ======================================================

# Usage:
# Rscript 04_spatial_marker_analysis.R <spaceranger_output_dir> <sample_id>

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 04_spatial_marker_analysis.R <spaceranger_output_dir> <sample_id>")
}

base_dir <- args[1]
sample_id <- args[2]

input_path <- file.path(
  base_dir,
  sample_id,
  "spatial_analysis",
  paste0(sample_id, "_spatial_clustered.rds")
)

if (!file.exists(input_path)) {
  stop("Clustered spatial object not found.")
}

output_dir <- file.path(base_dir, sample_id, "spatial_analysis", "markers")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("Loading clustered spatial object...\n")

spatial_obj <- readRDS(input_path)

# ------------------------------------------------------
# 1. Cluster Marker Identification
# ------------------------------------------------------

cat("Identifying spatial cluster markers...\n")

markers <- FindAllMarkers(
  spatial_obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(
  markers,
  file = file.path(output_dir, "spatial_cluster_markers.csv"),
  row.names = FALSE
)

# Top 10 per cluster
top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(
  top_markers,
  file = file.path(output_dir, "top10_spatial_markers_per_cluster.csv"),
  row.names = FALSE
)

# ------------------------------------------------------
# 2. Canonical Gene Visualization
# ------------------------------------------------------

cat("Generating canonical spatial feature plots...\n")

canonical_genes <- c("CD3D", "MS4A1", "LYZ", "COL1A1", "PECAM1")

pdf(file.path(output_dir, "canonical_spatial_features.pdf"))

for (gene in canonical_genes) {
  if (gene %in% rownames(spatial_obj)) {
    print(SpatialFeaturePlot(spatial_obj, features = gene))
  }
}

dev.off()

# ------------------------------------------------------
# 3. Region-Level Differential Scaffold
# ------------------------------------------------------

cat("Preparing region-level differential scaffold...\n")

# Example scaffold:
# Idents(spatial_obj) <- spatial_obj$region_label
# region_markers <- FindMarkers(spatial_obj, ident.1 = "RegionA", ident.2 = "RegionB")

saveRDS(
  spatial_obj,
  file = file.path(output_dir, "spatial_marker_ready_object.rds")
)

cat("Spatial marker analysis complete.\n")
