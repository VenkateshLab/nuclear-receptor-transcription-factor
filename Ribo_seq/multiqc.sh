#!/usr/bin/env bash
set -e

# Pre-trim QC
multiqc output/*/QC_check/ -o multiqc_pre_trim

# Post-trim QC
multiqc output/*/Fastp_fastqc/ -o multiqc_post_trim