#!/usr/bin/env bash
set -o errexit
set -o nounset

# Path to the data directory and QC script
DATA_DIR="data"
QC_SCRIPT="./QC.sh"

# Ensure the QC script is executable
chmod +x $QC_SCRIPT

# Read sample names from data.txt and process them
while read -r line; do
    # Extract the sample name (e.g., 1_S23 from 1_S23_R1.fastq.gz)
    SAMPLE_NAME=$(echo "$line" | sed 's/_R[12]\.fastq\.gz//')

    # Skip if this sample has already been processed (to avoid duplicates)
    if [[ ! -d "output/${SAMPLE_NAME}" ]]; then
        echo "Processing sample: $SAMPLE_NAME"
        $QC_SCRIPT "$SAMPLE_NAME" "$DATA_DIR"
    else
        echo "Sample $SAMPLE_NAME already processed. Skipping."
    fi
done < data.txt