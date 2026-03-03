# ======================================================
# Single-Cell RNA-seq: Cluster Annotation
# Omics Core Standardized Pipeline
# ======================================================

# Usage:
# Rscript 04_annotation.R <cellranger_output_dir> [annotation_mapping.csv]

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript 04_annotation.R <cellranger_output_dir> [annotation_mapping.csv]")
}

base_dir <- args[1]
annotation_file <- ifelse(length(args) >= 2, args[2], NA)

integrated_path <- file.path(base_dir, "integrated_analysis", "integrated_seurat_object.rds")

if (!file.exists(integrated_path)) {
  stop("Integrated Seurat object not found.")
}

obj <- readRDS(integrated_path)

output_dir <- file.path(base_dir, "integrated_analysis")
dir.create(output_dir, showWarnings = FALSE)

# ------------------------------------------------------
# Find Cluster Markers
# ------------------------------------------------------

markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(markers,
          file = file.path(output_dir, "all_cluster_markers.csv"),
          row.names = FALSE)

# ------------------------------------------------------
# Top Markers per Cluster
# ------------------------------------------------------

top_markers <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 10)

write.csv(top_markers,
          file = file.path(output_dir, "top10_markers_per_cluster.csv"),
          row.names = FALSE)

# ------------------------------------------------------
# Optional Manual Annotation
# ------------------------------------------------------

if (!is.na(annotation_file) && file.exists(annotation_file)) {

  cat("Applying manual annotation mapping...\n")

  annotation_map <- read.csv(annotation_file)

  # Expected format:
  # cluster,celltype

  obj$celltype <- plyr::mapvalues(
    x = as.character(Idents(obj)),
    from = as.character(annotation_map$cluster),
    to = annotation_map$celltype
  )

  Idents(obj) <- obj$celltype

  saveRDS(obj,
          file = file.path(output_dir, "annotated_seurat_object.rds"))

  pdf(file.path(output_dir, "UMAP_annotated.pdf"))
  DimPlot(obj, reduction = "umap", label = TRUE)
  dev.off()

} else {

  cat("No annotation mapping provided. Saving clustered object only.\n")

  saveRDS(obj,
          file = file.path(output_dir, "clustered_seurat_object.rds"))
}

cat("Annotation workflow complete.\n")
