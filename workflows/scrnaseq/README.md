# Single-Cell RNA-seq Workflow (10x Genomics)

This directory contains a standardized, modular pipeline for processing and analyzing 10x Genomics single-cell RNA-seq data.

The workflow is designed for core-facility use: configuration-driven, reproducible, and adaptable across projects.

---

## Workflow Overview

FASTQ → Cell Ranger → QC → Integration → Clustering → Annotation → Downstream Analysis

Each step is modular and can be executed independently.

---

## Pipeline Components

### 01_cellranger_processing.sh
Processes raw FASTQ files using 10x Genomics Cell Ranger.

**Input**
- Raw FASTQ files
- `config/sample_config.yaml`
- `metadata/sample_sheet.csv`

**Output**
- Cell Ranger output directories per sample
- Filtered feature-barcode matrices

---

### 02_qc_filtering.R
Performs initial Seurat object creation and quality control filtering.

**Tasks**
- Load filtered matrices
- Calculate mitochondrial percentage
- Apply standard filtering thresholds
- Save QC-filtered Seurat object
- Generate QC violin plots

**Output**
- `<sample_id>_qc_filtered.rds`
- QC plots

---

### 03_integration_clustering.R
Merges multiple samples and performs standard Seurat workflow.

**Tasks**
- Merge samples
- Normalize
- Identify variable features
- PCA
- UMAP
- Clustering

**Output**
- `integrated_seurat_object.rds`
- UMAP plots by sample and cluster

---

### 04_annotation.R
Identifies cluster markers and supports manual annotation.

**Tasks**
- Find cluster markers
- Export full marker table
- Export top markers per cluster
- Optional cluster-to-celltype mapping

**Output**
- `all_cluster_markers.csv`
- `top10_markers_per_cluster.csv`
- Annotated Seurat object
- Annotated UMAP

---

### 05_downstream_analysis.R
Advanced analysis module.

**Includes**
- Cluster-level differential testing
- Canonical marker visualization
- Receptor–ligand scaffold (CellChat-ready object)
- Trajectory scaffolds (Monocle3, Slingshot)

**Output**
- Marker tables
- Canonical marker plots
- Monocle3 CDS object
- Slingshot SCE object
- CellChat-ready Seurat object

---

## Required Software

- Cell Ranger
- R (>= 4.0)
- Seurat
- DESeq2
- monocle3
- slingshot
- SingleCellExperiment
- tidyverse

---

## Execution Order

```bash
# 1. Run Cell Ranger
bash workflows/scrnaseq/01_cellranger_processing.sh config/sample_config.yaml metadata/sample_sheet.csv

# 2. QC filtering (per sample)
Rscript workflows/scrnaseq/02_qc_filtering.R <cellranger_output_dir> Sample1

# 3. Integration & clustering
Rscript workflows/scrnaseq/03_integration_clustering.R <cellranger_output_dir> Sample1,Sample2

# 4. Annotation
Rscript workflows/scrnaseq/04_annotation.R <cellranger_output_dir>

# 5. Downstream analysis
Rscript workflows/scrnaseq/05_downstream_analysis.R <cellranger_output_dir>
