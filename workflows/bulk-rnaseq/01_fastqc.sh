#!/bin/bash

# ======================================================
# Bulk RNA-seq: FASTQ Quality Control
# Omics Core Standardized Pipeline
# ======================================================

# Usage:
# bash 01_fastqc.sh config/sample_config.yaml

set -e

CONFIG=$1

if [ -z "$CONFIG" ]; then
  echo "Error: Please provide a config file."
  echo "Usage: bash 01_fastqc.sh config/sample_config.yaml"
  exit 1
fi

# Extract parameters from config
RAW_DIR=$(grep raw_fastq_dir $CONFIG | awk '{print $2}')
OUT_DIR=$(grep qc_output_dir $CONFIG | awk '{print $2}')

if [ ! -d "$RAW_DIR" ]; then
  echo "Error: FASTQ directory not found."
  exit 1
fi

mkdir -p ${OUT_DIR}

echo "Running FastQC..."
for file in ${RAW_DIR}/*.fastq.gz
do
  fastqc $file -o ${OUT_DIR}
done

echo "FASTQ QC complete."
