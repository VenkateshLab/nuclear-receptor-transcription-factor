#!/usr/bin/env bash
set -o errexit
set -o nounset

# Reference genome file
REFERENCE="/home/manoj/Desktop/Cut_run/Reference/default/0f10d83b1050c08dd53189986f60970b92a315aa7a16a6f1.fa"

# Ensure the reference exists
if [ ! -f "$REFERENCE" ]; then
    echo "Reference genome not found: $REFERENCE"
    exit 1
fi

# Loop through each sample in data.txt
while read -r line; do
    # Extract the sample name (e.g., 1_S23 from 1_S23_R1.fastq.gz)
    SAMPLE_NAME=$(echo "$line" | sed 's/_R[12]\.fastq\.gz//' | uniq)

    # Skip if already processed
    if [[ -f "output/$SAMPLE_NAME/peaks/${SAMPLE_NAME}_peaks.narrowPeak" ]]; then
        echo "Sample $SAMPLE_NAME already processed. Skipping."
        continue
    fi

    echo "Processing sample: $SAMPLE_NAME"

    # Create output directories
    mkdir -p output/$SAMPLE_NAME/peaks

    # BWA alignment
    bwa mem -t 30 "$REFERENCE" \
    output/$SAMPLE_NAME/trimedreads_fastp/${SAMPLE_NAME}_1.fq.gz \
    output/$SAMPLE_NAME/trimedreads_fastp/${SAMPLE_NAME}_2.fq.gz | \
    samtools view -@ 30 -bS -o output/$SAMPLE_NAME/${SAMPLE_NAME}.bam -

    # Sort BAM file
    samtools sort -@ 30 output/$SAMPLE_NAME/${SAMPLE_NAME}.bam \
        -o output/$SAMPLE_NAME/${SAMPLE_NAME}_sort.bam

    # Query name sort and fixmate
    samtools sort -n -@ 30 -o output/$SAMPLE_NAME/${SAMPLE_NAME}_querysort.bam \
        output/$SAMPLE_NAME/${SAMPLE_NAME}_sort.bam
    samtools fixmate -@ 30 -m output/$SAMPLE_NAME/${SAMPLE_NAME}_querysort.bam \
        output/$SAMPLE_NAME/${SAMPLE_NAME}_fixmate.bam

    # Final sort and mark duplicates
    samtools sort -@ 30 -o output/$SAMPLE_NAME/${SAMPLE_NAME}_fixmate_sort.bam \
        output/$SAMPLE_NAME/${SAMPLE_NAME}_fixmate.bam
    samtools markdup -@ 30 -r output/$SAMPLE_NAME/${SAMPLE_NAME}_fixmate_sort.bam \
        output/$SAMPLE_NAME/${SAMPLE_NAME}_PCR.bam

    # Index the BAM file
    samtools index -@ 30 output/$SAMPLE_NAME/${SAMPLE_NAME}_PCR.bam

    # Generate BigWig file
    bamCoverage -p 30 -b output/$SAMPLE_NAME/${SAMPLE_NAME}_PCR.bam \
        -o output/$SAMPLE_NAME/${SAMPLE_NAME}.bw

    # Peak calling with MACS2
    macs2 callpeak -t output/$SAMPLE_NAME/${SAMPLE_NAME}_PCR.bam \
        -c output/$SAMPLE_NAME/${SAMPLE_NAME}_PCR.bam \
        --bdg --SPMR -n ${SAMPLE_NAME}_peaks \
        --outdir output/$SAMPLE_NAME/peaks

    echo "Finished processing sample: $SAMPLE_NAME"
done < data.txt
