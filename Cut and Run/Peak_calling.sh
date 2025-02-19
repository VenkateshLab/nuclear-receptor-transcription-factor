#!/usr/bin/env bash
set -o errexit
set -o nounset

while read -r line; do
    SAMPLE_NAME=$(echo "$line" | sed 's/_R[12]\.fastq\.gz//' | uniq)

    echo "Running peak calling for sample: $SAMPLE_NAME"

    mkdir -p output/$SAMPLE_NAME/peaks_new

    macs2 callpeak -t output/$SAMPLE_NAME/${SAMPLE_NAME}_PCR.bam \
        -f BAM -g mm -n ${SAMPLE_NAME} --nomodel --shift -100 --extsize 200 \
        --outdir output/$SAMPLE_NAME/peaks_new --keep-dup all
done < data.txt