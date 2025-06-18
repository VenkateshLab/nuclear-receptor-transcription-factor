#!/usr/bin/env bash
set -e

# ------------------------
# SET THESE VARIABLES!
# ------------------------
GTF_PATH="/media/manoj/Backup_3/Bulk_RNA/Ref/genes.gtf" 
THREADS=20
OUTPUT="featurecounts_output.txt"
# ------------------------

# Input BAMs
BAM_FILES=$(find output -type f -name "*_Aligned.sortedByCoord.out.bam" | tr '\n' ' ')

# Run featureCounts
featureCounts \
  -T ${THREADS} \
  -a "${GTF_PATH}" \
  -o "${OUTPUT}" \
  -g gene_name \
  -t exon \
  -p \
  -B \
  -C \
  ${BAM_FILES}

echo "✓ featureCounts completed: ${OUTPUT}"
