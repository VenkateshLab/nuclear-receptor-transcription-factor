#!/usr/bin/env bash
set -euo pipefail

# -------------
# Usage:
#   bash star_align.sh <sample> <star_index_dir> <threads>
#
# Example:
#   bash star_align.sh C41_S3_L004 /path/to/STARindex 16
# -------------

SAMPLE="$1"          # e.g. C41_S3_L004
STAR_INDEX="/media/manoj/Backup_3/Bulk_RNA/Ref/GenomeDir"      # absolute path to pre-built STAR genome index
THREADS="20"    # default to 8 if not provided

# Input FASTQ (trimmed by fastp)
R1="output/${SAMPLE}/trimedreads_fastp/${SAMPLE}_1.fq.gz"
R2="output/${SAMPLE}/trimedreads_fastp/${SAMPLE}_2.fq.gz"

# Output directories
ALIGN_DIR="output/${SAMPLE}/STAR"
mkdir -p "${ALIGN_DIR}"

# Run STAR
/media/manoj/Backup_3/Bulk_RNA/Ref/STAR/bin/Linux_x86_64/STAR \
  --runThreadN "${THREADS}" \
  --genomeDir "${STAR_INDEX}" \
  --readFilesIn "${R1}" "${R2}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${ALIGN_DIR}/${SAMPLE}_" \
  --outSAMtype BAM SortedByCoordinate \
  --outSAMattrRGline ID:${SAMPLE} SM:${SAMPLE} PL:ILLUMINA \
  --quantMode GeneCounts

# STAR creates <prefix>_Aligned.sortedByCoord.out.bam
# Optional index (if samtools available)
if command -v samtools &>/dev/null; then
  samtools index "${ALIGN_DIR}/${SAMPLE}_Aligned.sortedByCoord.out.bam"
fi

echo "✓ ${SAMPLE} alignment complete."
