# Bulk RNA-seq Workflow

This directory contains the standardized bulk RNA-seq processing pipeline used in the Omics Core infrastructure framework.

The workflow is modular and configuration-driven, allowing adaptation across projects while maintaining reproducibility.

---

## Workflow Overview

FASTQ → QC → Alignment → Quantification → Differential Expression

Scripts are executed sequentially using a shared configuration file.

---

## Pipeline Components

### 01_fastqc.sh
Performs quality control on raw FASTQ files using FastQC.

Input:
- Raw FASTQ files (paired-end)
- `config/sample_config.yaml`

Output:
- Per-sample FastQC reports

---

### 02_alignment.sh
Aligns reads to the reference genome using HISAT2 and generates sorted, indexed BAM files.

Input:
- Raw FASTQ files
- HISAT2 reference index

Output:
- Sorted BAM files
- BAM index files (.bai)

---

### 03_quantification.sh
Performs gene-level quantification using featureCounts.

Input:
- Sorted BAM files
- Annotation GTF file

Output:
- `gene_counts.txt`

---

### 04_differential_expression.R
Performs differential expression analysis using DESeq2.

Input:
- `gene_counts.txt`
- Metadata file (CSV with sample IDs and condition column)

Output:
- `differential_expression_results.csv`
- `normalized_counts.csv`
- `PCA_plot.pdf`

---

## Required Software

- FastQC
- HISAT2
- Samtools
- Subread (featureCounts)
- R (>= 4.0)
- DESeq2
- tidyverse

---

## Execution Order

```bash
bash workflows/bulk-rnaseq/01_fastqc.sh config/sample_config.yaml
bash workflows/bulk-rnaseq/02_alignment.sh config/sample_config.yaml
bash workflows/bulk-rnaseq/03_quantification.sh config/sample_config.yaml
Rscript workflows/bulk-rnaseq/04_differential_expression.R \
    count_output/gene_counts.txt \
    metadata/sample_metadata.csv
