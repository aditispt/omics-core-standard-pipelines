#!/bin/bash

# ======================================================
# Bulk RNA-seq: Gene-Level Quantification (featureCounts)
# Omics Core Standardized Pipeline
# ======================================================

# Usage:
# bash 03_quantification.sh config/sample_config.yaml

set -e

CONFIG=$1

if [ -z "$CONFIG" ]; then
  echo "Error: Please provide a config file."
  echo "Usage: bash 03_quantification.sh config/sample_config.yaml"
  exit 1
fi

# Extract parameters from config
ALIGN_DIR=$(grep alignment_output_dir $CONFIG | awk '{print $2}')
COUNT_DIR=$(grep count_output_dir $CONFIG | awk '{print $2}')
GTF=$(grep annotation_gtf $CONFIG | awk '{print $2}')
THREADS=$(grep threads $CONFIG | awk '{print $2}')

mkdir -p ${COUNT_DIR}

echo "Running featureCounts..."

featureCounts \
  -T ${THREADS} \
  -a ${GTF} \
  -o ${COUNT_DIR}/gene_counts.txt \
  ${ALIGN_DIR}/*.sorted.bam

echo "Quantification complete."
