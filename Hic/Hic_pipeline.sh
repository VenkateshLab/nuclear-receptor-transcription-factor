#!/bin/sh

# Automated Hi-C Pipeline using Homer
# This script reads sample information from data.txt and processes all samples

# Set reference genome path
REFERENCE="Reference/Bowtie/mm10/mm10"

# Set number of CPU cores to use
CPU=20

# Check if data.txt exists
if [ ! -f "data.txt" ]; then
    echo "Error: data.txt file not found!"
    exit 1
fi

# Function to process a single Hi-C sample
process_sample() {
    local r1="$1"
    local r2="$2"
    sample_name=$(echo "$r1" | sed 's/_R[12].*//g' | sed 's/_S[0-9]*//g')
    
    echo "Processing sample: $sample_name"
    
    # Create output directory
    mkdir -p "$sample_name"
    
    # Step 1: Trim reads
    echo "Step 1: Trimming reads..."
    homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 "$r1"
    homerTools trim -3 GATC -mis 0 -matchStart 20 -min 20 "$r2"
    
    r1_trimmed="${r1}.trimmed"
    r2_trimmed="${r2}.trimmed"
    
    # Step 2: Align reads to reference genome
    echo "Step 2: Aligning reads to reference genome..."
    bowtie2 -p $CPU -x $REFERENCE -U "$r1_trimmed" > "${sample_name}_R1_mm10.sam"
    bowtie2 -p $CPU -x $REFERENCE -U "$r2_trimmed" > "${sample_name}_R2_mm10.sam"
    
    # Step 3: Create tag directory
    echo "Step 3: Creating tag directory..."
    makeTagDirectory "${sample_name}_tagdir" "${sample_name}_R1_mm10.sam,${sample_name}_R2_mm10.sam" -tbp 1
    
    # Step 4: Run Hi-C PCA
    echo "Step 4: Running Hi-C PCA..."
    runHiCpca.pl auto "${sample_name}_tagdir/" -res 25000 -window 50000 -genome mm10 -cpu $CPU
    
    # Step 5: Analyze Hi-C data
    echo "Step 5: Analyzing Hi-C data..."
    analyzeHiC "${sample_name}_tagdir/" -res 5000 -window 15000 -nomatrix -compactionStats auto -cpu $CPU
    
    # Step 6: Find TADs and loops
    echo "Step 6: Finding TADs and loops..."
    findTADsAndLoops.pl find "${sample_name}_tagdir/" -cpu 10 -res 3000 -window 15000 -genome mm10
    
    echo "Processing completed for $sample_name"
    echo "=========================================="
}

# Main execution
echo "Starting Hi-C pipeline..."

# Create a temporary file to store sample pairs
temp_file=$(mktemp)

# Parse data.txt file and find pairs
grep "_R1" data.txt > "$temp_file.r1"
grep "_R2" data.txt > "$temp_file.r2"

# Process each R1 file and find its matching R2
cat "$temp_file.r1" | while read r1_file; do
    # Extract base name (remove _R1 part)
    base=$(echo "$r1_file" | sed 's/_R1.*//g')
    
    # Find matching R2 file
    r2_file=$(grep "${base}_R2" "$temp_file.r2" | head -1)
    
    # Check if we have both files and they exist
    if [ -n "$r1_file" ] && [ -n "$r2_file" ] && [ -f "$r1_file" ] && [ -f "$r2_file" ]; then
        echo "Found sample pair: $r1_file and $r2_file"
        process_sample "$r1_file" "$r2_file"
    else
        echo "Warning: Could not find complete pair for sample $base. Skipping."
    fi
done

# Clean up temporary files
rm -f "$temp_file" "$temp_file.r1" "$temp_file.r2"

echo "Hi-C pipeline completed for all samples in data.txt!"
