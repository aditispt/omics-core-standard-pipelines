# ======================================================
# Single-Cell RNA-seq: Downstream Analysis Module
# Omics Core Standardized Pipeline
# ======================================================

# Usage:
# Rscript 05_downstream_analysis.R <cellranger_output_dir>

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop("Usage: Rscript 05_downstream_analysis.R <cellranger_output_dir>")
}

base_dir <- args[1]

integrated_path <- file.path(base_dir, "integrated_analysis", "annotated_seurat_object.rds")

if (!file.exists(integrated_path)) {
  integrated_path <- file.path(base_dir, "integrated_analysis", "clustered_seurat_object.rds")
}

if (!file.exists(integrated_path)) {
  stop("No integrated Seurat object found.")
}

obj <- readRDS(integrated_path)

output_dir <- file.path(base_dir, "integrated_analysis", "downstream")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------
# 1. Cluster-level Differential Expression
# ------------------------------------------------------

cat("Running cluster-level differential testing...\n")

cluster_markers <- FindAllMarkers(
  obj,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

write.csv(cluster_markers,
          file = file.path(output_dir, "cluster_level_markers.csv"),
          row.names = FALSE)

# ------------------------------------------------------
# 2. Canonical Marker Testing
# ------------------------------------------------------

cat("Testing canonical markers...\n")

canonical_markers <- c("CD3D", "CD4", "CD8A", "MS4A1", "LYZ", "NKG7")

pdf(file.path(output_dir, "canonical_marker_expression.pdf"))
FeaturePlot(obj, features = canonical_markers)
dev.off()

# ------------------------------------------------------
# 3. Receptor-Ligand Scaffold (CellChat-ready)
# ------------------------------------------------------

cat("Preparing object for receptor-ligand analysis...\n")

# This section prepares the object; users can extend with CellChat or similar
DefaultAssay(obj) <- "RNA"
saveRDS(obj,
        file = file.path(output_dir, "cellchat_ready_object.rds"))

# ------------------------------------------------------
# 4. Trajectory Analysis (Monocle3 Scaffold)
# ------------------------------------------------------

cat("Preparing Monocle3 trajectory scaffold...\n")

suppressPackageStartupMessages(library(monocle3))

cds <- as.cell_data_set(obj)

cds <- cluster_cells(cds)
cds <- learn_graph(cds)

saveRDS(cds,
        file = file.path(output_dir, "monocle3_cds_object.rds"))

# ------------------------------------------------------
# 5. Trajectory Analysis (Slingshot Scaffold)
# ------------------------------------------------------

cat("Preparing Slingshot trajectory scaffold...\n")

suppressPackageStartupMessages(library(slingshot))

sce <- as.SingleCellExperiment(obj)

sce <- slingshot(sce, clusterLabels = Idents(obj), reducedDim = "PCA")

saveRDS(sce,
        file = file.path(output_dir, "slingshot_sce_object.rds"))

cat("Downstream analysis module complete.\n")
