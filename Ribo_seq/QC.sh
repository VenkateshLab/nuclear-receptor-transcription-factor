#!/usr/bin/env bash
set -o errexit
set -o nounset

# $1 = sample name (e.g., C41_S3_L004)
# $2 = working directory (not directly used, but accepted for compatibility)

# Create output directories
mkdir -p output/${1}/QC_check
mkdir -p output/${1}/trimedreads_fastp
mkdir -p output/${1}/Fastp_fastqc

# Step 1: FastQC on raw reads
fastqc -t 30 ${1}_R1_001.fastq.gz -o output/${1}/QC_check
fastqc -t 30 ${1}_R2_001.fastq.gz -o output/${1}/QC_check

# Step 2: Trimming with fastp
fastp \
  -i ${1}_R1_001.fastq.gz \
  -I ${1}_R2_001.fastq.gz \
  -o output/${1}/trimedreads_fastp/${1}_1.fq.gz \
  -O output/${1}/trimedreads_fastp/${1}_2.fq.gz \
  -j output/${1}/trimedreads_fastp/${1}.json \
  -h output/${1}/trimedreads_fastp/${1}.html \
  -w 16

# Step 3: FastQC on trimmed reads
fastqc -t 30 output/${1}/trimedreads_fastp/${1}_1.fq.gz -o output/${1}/Fastp_fastqc
fastqc -t 30 output/${1}/trimedreads_fastp/${1}_2.fq.gz -o output/${1}/Fastp_fastqc

