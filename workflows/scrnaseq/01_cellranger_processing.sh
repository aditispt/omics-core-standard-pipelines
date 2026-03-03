#!/bin/bash

# ======================================================
# Single-Cell RNA-seq: 10x Genomics Processing (Cell Ranger)
# Omics Core Standardized Pipeline
# ======================================================

# Usage:
# bash 01_cellranger_processing.sh config/sample_config.yaml metadata/sample_sheet.csv

set -e

CONFIG=$1
SAMPLE_SHEET=$2

if [ -z "$CONFIG" ] || [ -z "$SAMPLE_SHEET" ]; then
  echo "Usage: bash 01_cellranger_processing.sh config/sample_config.yaml metadata/sample_sheet.csv"
  exit 1
fi

# ------------------------------------------------------
# Extract parameters from config
# ------------------------------------------------------

RAW_FASTQ_DIR=$(grep raw_fastq_dir $CONFIG | awk '{print $2}')
CELLRANGER_OUTPUT_DIR=$(grep cellranger_output_dir $CONFIG | awk '{print $2}')
REFERENCE_TRANSCRIPTOME=$(grep cellranger_reference $CONFIG | awk '{print $2}')
THREADS=$(grep threads $CONFIG | awk '{print $2}')
MEMORY=$(grep memory_gb $CONFIG | awk '{print $2}')

mkdir -p ${CELLRANGER_OUTPUT_DIR}

echo "Starting Cell Ranger processing..."

# ------------------------------------------------------
# Loop through sample sheet
# Sample sheet format:
# sample_id,fastq_subdir
# ------------------------------------------------------

while IFS=',' read -r SAMPLE_ID FASTQ_SUBDIR
do
  if [ "$SAMPLE_ID" == "sample_id" ]; then
    continue
  fi

  echo "Processing sample: ${SAMPLE_ID}"

  cellranger count \
    --id=${SAMPLE_ID} \
    --transcriptome=${REFERENCE_TRANSCRIPTOME} \
    --fastqs=${RAW_FASTQ_DIR}/${FASTQ_SUBDIR} \
    --sample=${SAMPLE_ID} \
    --localcores=${THREADS} \
    --localmem=${MEMORY}

  mv ${SAMPLE_ID} ${CELLRANGER_OUTPUT_DIR}/

done < ${SAMPLE_SHEET}

echo "Cell Ranger processing complete."
