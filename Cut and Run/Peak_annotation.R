#!/usr/bin/env Rscript

# Load necessary libraries
library(ChIPseeker)
library(GenomicRanges)
library(org.Mm.eg.db)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(data.table)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if correct number of arguments are provided
if (length(args) != 2) {
  stop("Usage: Rscript myscript.R inputfile.xls outputfile.csv")
}

# Assign arguments to variables
input_file <- args[1]
output_file <- args[2]

# Read the file, manually filtering out lines that start with '#'
lines <- readLines(input_file)
data_lines <- lines[!grepl("^#", lines)]  # Remove lines starting with #

# Write the filtered data to a temporary file
temp_file <- tempfile()
writeLines(data_lines, temp_file)

# Load the cleaned data into R
peak_data <- fread(temp_file, sep="\t", header=TRUE, stringsAsFactors=FALSE)

# Convert to GenomicRanges object
peaks <- makeGRangesFromDataFrame(df = peak_data, 
                                  keep.extra.columns = TRUE,
                                  seqnames.field = "chr", 
                                  start.field = "start", 
                                  end.field = "end")

# Annotate peaks
peakAnno <- annotatePeak(peaks, 
                         TxDb = TxDb.Mmusculus.UCSC.mm10.knownGene, 
                         annoDb = "org.Mm.eg.db")

# Convert annotation to data frame
peakAnno_df <- as.data.frame(peakAnno)

# Save to CSV
write.csv(peakAnno_df, file = output_file, row.names = FALSE)

cat("Annotation completed. Output saved to:", output_file, "\n")
