# Spatial Transcriptomics Workflow (10x Visium)

This directory contains a standardized, modular workflow for processing and analyzing 10x Genomics Visium spatial transcriptomics data.

The pipeline is designed for core-facility deployment: configuration-driven, reproducible, and adaptable across tissue types and biological questions.

---

## Workflow Overview

FASTQ → Space Ranger → QC & Normalization → Clustering → Spatial Mapping → Marker & Region Analysis

Each step is modular and can be executed independently.

---

## Pipeline Components

### 01_spaceranger_processing.sh
Processes raw FASTQ files using 10x Genomics Space Ranger.

**Input**
- Raw FASTQ files
- Spatial configuration file
- `metadata/sample_sheet.csv` (sample_id, fastq_subdir, image_path)

**Output**
- Space Ranger output directory per sample
- Filtered feature-barcode matrix
- Spatial tissue position files

---

### 02_qc_normalization.R
Performs Seurat spatial object creation and quality control.

**Tasks**
- Load Space Ranger output
- Calculate mitochondrial percentage
- Apply conservative filtering thresholds
- Normalize using SCTransform
- Generate spatial QC plots

**Output**
- `<sample_id>_spatial_normalized.rds`
- Spatial QC plots

---

### 03_clustering_spatial_mapping.R
Performs dimensionality reduction and clustering.

**Tasks**
- PCA
- Neighbor graph construction
- UMAP embedding
- Clustering
- Spatial cluster visualization

**Output**
- Clustered spatial Seurat object
- UMAP plot
- Spatial cluster map

---

### 04_spatial_marker_analysis.R
Performs spatial cluster marker analysis and region-level scaffolding.

**Tasks**
- Identify cluster markers
- Export full marker table
- Export top markers per cluster
- Visualize canonical genes in spatial context
- Provide scaffold for region-based differential testing

**Output**
- `spatial_cluster_markers.csv`
- `top10_spatial_markers_per_cluster.csv`
- Canonical spatial feature plots
- Spatial marker-ready object

---

## Required Software

- Space Ranger
- R (>= 4.0)
- Seurat
- tidyverse
- ggplot2

---

## Execution Order

```bash
# 1. Run Space Ranger
bash workflows/spatial-transcriptomics/01_spaceranger_processing.sh \
    config/spatial_config.yaml \
    metadata/sample_sheet.csv

# 2. QC & normalization
Rscript workflows/spatial-transcriptomics/02_qc_normalization.R \
    <spaceranger_output_dir> Sample1

# 3. Clustering & spatial mapping
Rscript workflows/spatial-transcriptomics/03_clustering_spatial_mapping.R \
    <spaceranger_output_dir> Sample1

# 4. Marker analysis
Rscript workflows/spatial-transcriptomics/04_spatial_marker_analysis.R \
    <spaceranger_output_dir> Sample1
