#!/bin/bash

# Create an output directory for trimmed reads
mkdir -p trimmed_reads

# Loop through R1 FASTQ files and find corresponding R2 files
for R1 in *_R1_001.fastq.gz; do
    # Derive R2 filename
    R2=${R1/_R1_001.fastq.gz/_R2_001.fastq.gz}

    # Check if R2 file exists
    if [[ -f "$R2" ]]; then
        # Define output filenames
        SAMPLE_NAME=$(basename "$R1" _R1_001.fastq.gz)
        OUTPUT_R1="trimmed_reads/${SAMPLE_NAME}_trimmed_R1.fastq.gz"
        OUTPUT_R2="trimmed_reads/${SAMPLE_NAME}_trimmed_R2.fastq.gz"
        REPORT_HTML="trimmed_reads/${SAMPLE_NAME}_fastp_report.html"
        REPORT_JSON="trimmed_reads/${SAMPLE_NAME}_fastp_report.json"

        # Run fastp
        fastp \
            --in1 "$R1" --in2 "$R2" \
            --out1 "$OUTPUT_R1" --out2 "$OUTPUT_R2" \
            --qualified_quality_phred 30 \
            --detect_adapter_for_pe \
            --length_required 20 \
            --html "$REPORT_HTML" \
            --json "$REPORT_JSON" \
            --thread 10

        echo "Processed: $SAMPLE_NAME"
    else
        echo "Skipping $R1 (No matching R2 file found)"
    fi
done

echo "Fastp trimming complete. Trimmed reads are in 'trimmed_reads/'"