#!/usr/bin/env bash
set -o errexit
set -o nounset
set -o pipefail

# Reference genome file
REFERENCE="/home/manoj/Desktop/Cut_run/Reference/default/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1.fa"

# Ensure the reference exists
if [ ! -f "$REFERENCE" ]; then
    echo "Reference genome not found: $REFERENCE"
    exit 1
fi

# Loop through each sample in data.txt
while read -r line; do
    SAMPLE_NAME=$(echo "$line" | sed 's/_R[12]\.fastq\.gz//' | uniq)

    if [ -f "output/$SAMPLE_NAME/${SAMPLE_NAME}_PCR.bam" ]; then
        echo "Sample $SAMPLE_NAME already processed. Skipping."
        continue
    fi

    echo "Processing sample: $SAMPLE_NAME"
    
    # BWA alignment
    bwa mem -t 30 "$REFERENCE" \
        output/$SAMPLE_NAME/trimedreads_fastp/${SAMPLE_NAME}_1.fq.gz \
        output/$SAMPLE_NAME/trimedreads_fastp/${SAMPLE_NAME}_2.fq.gz | \
    samtools view -@ 30 -bS -o output/$SAMPLE_NAME/${SAMPLE_NAME}.bam - || {
        echo "Error in alignment for sample $SAMPLE_NAME"
        continue
    }

    # Sort BAM file
    samtools sort -@ 30 output/$SAMPLE_NAME/${SAMPLE_NAME}.bam \
        -o output/$SAMPLE_NAME/${SAMPLE_NAME}_sort.bam || {
        echo "Error in sorting for sample $SAMPLE_NAME"
        continue
    }

    # Query name sort and fixmate
    samtools sort -n -@ 30 -o output/$SAMPLE_NAME/${SAMPLE_NAME}_querysort.bam \
        output/$SAMPLE_NAME/${SAMPLE_NAME}_sort.bam
    samtools fixmate -@ 30 -m output/$SAMPLE_NAME/${SAMPLE_NAME}_querysort.bam \
        output/$SAMPLE_NAME/${SAMPLE_NAME}_fixmate.bam || {
        echo "Error in fixmate for sample $SAMPLE_NAME"
        continue
    }

    # Final sort and mark duplicates
    samtools sort -@ 30 -o output/$SAMPLE_NAME/${SAMPLE_NAME}_fixmate_sort.bam \
        output/$SAMPLE_NAME/${SAMPLE_NAME}_fixmate.bam
    samtools markdup -@ 30 -r output/$SAMPLE_NAME/${SAMPLE_NAME}_fixmate_sort.bam \
        output/$SAMPLE_NAME/${SAMPLE_NAME}_PCR.bam || {
        echo "Error in marking duplicates for sample $SAMPLE_NAME"
        continue
    }

    # Index the BAM file
    samtools index -@ 30 output/$SAMPLE_NAME/${SAMPLE_NAME}_PCR.bam || {
        echo "Error in indexing for sample $SAMPLE_NAME"
        continue
    }

    # Cleanup intermediate files (optional)
    rm -f output/$SAMPLE_NAME/${SAMPLE_NAME}.bam \
        output/$SAMPLE_NAME/${SAMPLE_NAME}_sort.bam \
        output/$SAMPLE_NAME/${SAMPLE_NAME}_querysort.bam \
        output/$SAMPLE_NAME/${SAMPLE_NAME}_fixmate.bam \
        output/$SAMPLE_NAME/${SAMPLE_NAME}_fixmate_sort.bam

    echo "Finished processing sample: $SAMPLE_NAME"
done < data.txt
