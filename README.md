# omics-core-standard-pipelines

Standardized, forkable multi-omics workflow templates for bulk RNA-seq, single-cell RNA-seq, CITE-seq, and spatial transcriptomics.

This repository provides modular directory structures, execution templates, configuration-driven workflows, and reproducibility scaffolding designed for bioinformatics cores and translational research programs.

---

## Supported Modalities

### Bulk RNA-seq
FASTQ → QC → Alignment → Quantification → Differential Expression → Visualization  

### Single-Cell RNA-seq
Preprocessing → QC → Normalization → Clustering → Annotation → Differential Expression  

### CITE-seq
RNA + ADT integration → Weighted Nearest Neighbor analysis → Cell-type scoring  

### Spatial Transcriptomics
Alignment → Spatial clustering → Feature visualization → Region-level analysis  

---

## Repository Structure

config/  
User-editable configuration files  

workflows/  
Pipeline modules organized by modality  

environments/  
Reproducible software environment definitions  

templates/  
Metadata and sample sheet templates  

---

Each workflow stage is modular and can be executed independently or integrated into workflow managers such as Snakemake or Nextflow.
