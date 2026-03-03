#!/bin/bash

# ======================================================
# Spatial Transcriptomics: 10x Visium Processing
# Omics Core Standardized Pipeline
# ======================================================

# Usage:
# bash 01_spaceranger_processing.sh config/spatial_config.yaml metadata/sample_sheet.csv

set -e

CONFIG=$1
SAMPLE_SHEET=$2

if [ -z "$CONFIG" ] || [ -z "$SAMPLE_SHEET" ]; then
  echo "Usage: bash 01_spaceranger_processing.sh config/spatial_config.yaml metadata/sample_sheet.csv"
  exit 1
fi

RAW_FASTQ_DIR=$(grep raw_fastq_dir $CONFIG | awk '{print $2}')
SPACERANGER_OUTPUT_DIR=$(grep spaceranger_output_dir $CONFIG | awk '{print $2}')
REFERENCE_TRANSCRIPTOME=$(grep spaceranger_reference $CONFIG | awk '{print $2}')
THREADS=$(grep threads $CONFIG | awk '{print $2}')
MEMORY=$(grep memory_gb $CONFIG | awk '{print $2}')

mkdir -p ${SPACERANGER_OUTPUT_DIR}

echo "Starting Space Ranger processing..."

while IFS=',' read -r SAMPLE_ID FASTQ_SUBDIR IMAGE_PATH
do
  if [ "$SAMPLE_ID" == "sample_id" ]; then
    continue
  fi

  echo "Processing spatial sample: ${SAMPLE_ID}"

  spaceranger count \
    --id=${SAMPLE_ID} \
    --transcriptome=${REFERENCE_TRANSCRIPTOME} \
    --fastqs=${RAW_FASTQ_DIR}/${FASTQ_SUBDIR} \
    --sample=${SAMPLE_ID} \
    --image=${IMAGE_PATH} \
    --localcores=${THREADS} \
    --localmem=${MEMORY}

  mv ${SAMPLE_ID} ${SPACERANGER_OUTPUT_DIR}/

done < ${SAMPLE_SHEET}

echo "Space Ranger processing complete."
