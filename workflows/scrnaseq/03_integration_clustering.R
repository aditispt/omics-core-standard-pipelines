# ======================================================
# Single-Cell RNA-seq: Integration & Clustering (Seurat)
# Omics Core Standardized Pipeline
# ======================================================

# Usage:
# Rscript 03_integration_clustering.R <cellranger_output_dir> <sample1,sample2,...>

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Rscript 03_integration_clustering.R <cellranger_output_dir> <sample1,sample2,...>")
}

base_dir <- args[1]
sample_list <- strsplit(args[2], ",")[[1]]

seurat_list <- list()

# ------------------------------------------------------
# Load QC-filtered objects
# ------------------------------------------------------

for (sample_id in sample_list) {

  input_path <- file.path(
    base_dir,
    sample_id,
    "seurat_objects",
    paste0(sample_id, "_qc_filtered.rds")
  )

  if (!file.exists(input_path)) {
    stop(paste("QC-filtered file not found for sample:", sample_id))
  }

  obj <- readRDS(input_path)
  obj$orig.ident <- sample_id
  seurat_list[[sample_id]] <- obj
}

cat("Loaded", length(seurat_list), "samples\n")

# ------------------------------------------------------
# Merge Samples
# ------------------------------------------------------

combined <- merge(seurat_list[[1]],
                  y = seurat_list[-1],
                  add.cell.ids = sample_list)

# ------------------------------------------------------
# Standard Seurat Workflow
# ------------------------------------------------------

combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined, npcs = 30)
combined <- RunUMAP(combined, dims = 1:20)
combined <- FindNeighbors(combined, dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)

# ------------------------------------------------------
# Save Integrated Object
# ------------------------------------------------------

output_dir <- file.path(base_dir, "integrated_analysis")
dir.create(output_dir, showWarnings = FALSE)

saveRDS(combined,
        file = file.path(output_dir, "integrated_seurat_object.rds"))

# ------------------------------------------------------
# UMAP Plot
# ------------------------------------------------------

pdf(file.path(output_dir, "UMAP_by_sample.pdf"))
DimPlot(combined, reduction = "umap", group.by = "orig.ident")
dev.off()

pdf(file.path(output_dir, "UMAP_by_cluster.pdf"))
DimPlot(combined, reduction = "umap", label = TRUE)
dev.off()

cat("Integration and clustering complete.\n")
