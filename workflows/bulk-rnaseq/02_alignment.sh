#!/bin/bash

# ======================================================
# Bulk RNA-seq: Alignment (HISAT2)
# Omics Core Standardized Pipeline
# ======================================================

# Usage:
# bash 02_alignment.sh config/sample_config.yaml

set -e

CONFIG=$1

if [ -z "$CONFIG" ]; then
  echo "Error: Please provide a config file."
  echo "Usage: bash 02_alignment.sh config/sample_config.yaml"
  exit 1
fi

# Extract parameters from config
RAW_DIR=$(grep raw_fastq_dir $CONFIG | awk '{print $2}')
ALIGN_DIR=$(grep alignment_output_dir $CONFIG | awk '{print $2}')
REFERENCE_INDEX=$(grep reference_index $CONFIG | awk '{print $2}')
THREADS=$(grep threads $CONFIG | awk '{print $2}')

mkdir -p ${ALIGN_DIR}

echo "Starting alignment with HISAT2..."

for R1 in ${RAW_DIR}/*_R1*.fastq.gz
do
  R2=${R1/_R1/_R2}
  SAMPLE=$(basename ${R1} | cut -d "_" -f1)

  hisat2 \
    -x ${REFERENCE_INDEX} \
    -1 ${R1} \
    -2 ${R2} \
    -p ${THREADS} \
  | samtools sort -@ ${THREADS} -o ${ALIGN_DIR}/${SAMPLE}.sorted.bam

  samtools index ${ALIGN_DIR}/${SAMPLE}.sorted.bam

done

echo "Alignment complete."
