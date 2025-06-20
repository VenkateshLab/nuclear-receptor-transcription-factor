# Multi-Omics Analysis Pipeline

A comprehensive bioinformatics pipeline for analyzing multiple types of sequencing data including ATAC-seq, single-cell RNA-seq, Hi-C, bulk RNA-seq, and quality control analysis.

## **Overview**

This pipeline provides end-to-end analysis workflows for:

- **🧬 ATAC-seq**: Chromatin accessibility analysis using TOBIAS
- **🔬 Single-cell RNA-seq**: Cell Ranger and Seurat-based analysis
- **🌐 Hi-C**: 3D chromatin structure analysis using Homer
- **📊 Bulk RNA-seq**: STAR alignment and feature counting
- **✅ Quality Control**: FastQC and FastP preprocessing

---

## **Quick Start**

### **1. Clone Repository**
```bash
git clone <your-repository-url>
cd <your-repository-name>
chmod +x *.sh  # Make all shell scripts executable
```

### **2. Prepare Input Data**
Create a `data.txt` file listing your sample files:
```
sample_1_R1.fastq.gz
sample_1_R2.fastq.gz
sample_2_R1.fastq.gz
sample_2_R2.fastq.gz
```

### **3. Run Analysis**
```bash
# Quality Control
./QC.sh sample_1 ./

# ATAC-seq Analysis
./alignment.sh
./Peak_calling.sh

# Single-cell RNA-seq
./CellRanger.sh sample_1

# Hi-C Analysis  
./Hic_pipeline.sh

# Bulk RNA-seq
./star_align.sh sample_1
./feature_counts.sh
```

---

## **Required Tools and Dependencies**

### **Core Bioinformatics Tools**

| Tool | Purpose | Installation |
|------|---------|-------------|
| **FastQC** | Quality assessment | `conda install -c bioconda fastqc` |
| **FastP** | Read trimming | `conda install -c bioconda fastp` |
| **BWA** | DNA alignment | `conda install -c bioconda bwa` |
| **SAMtools** | BAM file processing | `conda install -c bioconda samtools` |
| **MACS2** | Peak calling | `conda install -c bioconda macs2` |
| **STAR** | RNA-seq alignment | `conda install -c bioconda star` |
| **featureCounts** | Gene quantification | `conda install -c bioconda subread` |
| **MultiQC** | Report aggregation | `conda install -c bioconda multiqc` |

### **Specialized Analysis Tools**

| Tool | Purpose | Installation |
|------|---------|-------------|
| **TOBIAS** | ATAC-seq footprinting | `conda install -c bioconda tobias` |
| **Homer** | Hi-C analysis | [Manual installation](http://homer.ucsd.edu/homer/download.html) |
| **Cell Ranger** | Single-cell processing | [10x Genomics download](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest) |
| **Docker** | Containerization | [Docker installation](https://docs.docker.com/get-docker/) |

### **R Libraries**

```r
# Essential R packages
install.packages(c("Seurat", "ggplot2", "dplyr", "plyr", "RColorBrewer"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "ChIPseeker", 
    "GenomicRanges", 
    "org.Mm.eg.db", 
    "TxDb.Mmusculus.UCSC.mm10.knownGene",
    "scCustomize",
    "ShinyCell"
))
```

---

## **Installation Guide**

### **Option 1: Conda Environment (Recommended)**

```bash
# Create conda environment
conda create -n multiomics python=3.8
conda activate multiomics

# Install bioinformatics tools
conda install -c bioconda -c conda-forge \
    fastqc fastp bwa samtools macs2 star subread \
    bowtie2 tobias snakemake multiqc

# Install Python packages
pip install tobias

# Install R
conda install -c conda-forge r-base r-essentials
```

### **Option 2: Docker Setup**

```bash
# Pull Cell Ranger Docker image
docker pull litd/docker-cellranger

# Verify Docker installation
docker run --rm litd/docker-cellranger cellranger --version
```

### **Option 3: Manual Installation**

```bash
# Install Homer for Hi-C analysis
wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install
perl configureHomer.pl -install mm10

# Download Cell Ranger
# Visit: https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest
```

---

## **Reference Genome Setup**

### **Download Required Genomes**

```bash
# Create reference directory
mkdir -p Reference

# Download mm10 genome for ATAC-seq/Hi-C
cd Reference
wget http://hgdownload.cse.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip mm10.fa.gz

# Index genome for BWA
bwa index mm10.fa

# Download Cell Ranger reference (for single-cell)
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
tar -zxvf refdata-gex-mm10-2020-A.tar.gz

# Download STAR genome index (for bulk RNA-seq)
# Build STAR index or download pre-built index
```

### **Required Reference Files**

- **FASTA**: `mm10.fa` (mouse genome)
- **GTF**: `genes.gtf` (gene annotations)
- **Blacklist**: `mm10_blacklist.bed` (for ATAC-seq)
- **Cell Ranger Reference**: Pre-built transcriptome reference

---

## **Pipeline Workflows**

### **1. Quality Control Workflow**

```bash
# Run QC on all samples
while read -r line; do
    SAMPLE_NAME=$(echo "$line" | sed 's/_R[12]\.fastq\.gz//' | uniq)
    ./QC.sh "$SAMPLE_NAME" ./
done < data.txt

# Generate MultiQC reports
./multiqc.sh
```

**Output:**
- `output/{sample}/QC_check/`: FastQC reports
- `output/{sample}/trimedreads_fastp/`: Trimmed reads
- `multiqc_pre_trim/`: Pre-trimming QC report
- `multiqc_post_trim/`: Post-trimming QC report

### **2. ATAC-seq Analysis Workflow**

```bash
# Step 1: Alignment
./alignment.sh

# Step 2: Peak calling
./Peak_calling.sh

# Step 3: Peak annotation
Rscript Peak_annotation.R input_peaks.xls output_annotated.csv

# Step 4: TOBIAS footprinting
./Tobias.sh
```

**Output:**
- `output/{sample}/{sample}_PCR.bam`: Processed alignments
- `output/{sample}/peaks_new/`: MACS2 peak calls
- `P0_P4_P7/`: TOBIAS footprinting results

### **3. Single-cell RNA-seq Workflow**

```bash
# Process with Cell Ranger
./CellRanger.sh sample_1

# Analyze with Seurat (R)
Rscript -e "rmarkdown::render('sc_RNA.Rmd')"
```

**Output:**
- `sample_1/outs/`: Cell Ranger outputs
- `NR3F1_F6_data.rds`: Seurat object
- Various plots and analysis results

### **4. Hi-C Analysis Workflow**

```bash
# Run Hi-C pipeline
./Hic_pipeline.sh
```

**Output:**
- `{sample}_tagdir/`: Homer tag directories
- TAD and loop calling results
- Hi-C interaction matrices

### **5. Bulk RNA-seq Workflow**

```bash
# Step 1: STAR alignment
./star_align.sh sample_1 /path/to/STARindex 20

# Step 2: Feature counting
./feature_counts.sh

# Step 3: Generate reports
./multiqc.sh
```

**Output:**
- `output/{sample}/STAR/`: STAR alignment results
- `featurecounts_output.txt`: Gene count matrix

---

## **Configuration Files**

### **TOBIAS Configuration (`config.yaml`)**

```yaml
data:
  P0: [/path/to/P0.bam]
  P4: [/path/to/P4.bam]
  P7: [/path/to/P7.bam]

run_info:
  organism: mouse
  fasta: /path/to/mm10.fa
  blacklist: /path/to/mm10_blacklist.bed
  gtf: /path/to/mm10_v2.gtf
  motifs: /path/to/Jaspar/*
  output: P0_P4_P7
```

### **Cell Ranger Parameters**

- **Reference Genome**: `Reference_genome_transgene/Transgene_Reference/`
- **Local Memory**: 200 GB
- **Local Cores**: 40

---

## **Input File Requirements**

### **File Structure**

```
project/
├── data.txt                     # Sample list
├── Raw_data/                    # Input FASTQ files
│   ├── sample_1_R1.fastq.gz
│   ├── sample_1_R2.fastq.gz
│   └── ...
├── Reference/                   # Reference genomes
├── config.yaml                  # TOBIAS configuration
└── scripts/                     # Analysis scripts
```

### **Sample Naming Convention**

FASTQ files should follow the pattern:
- `{sample}_R1.fastq.gz` (forward reads)
- `{sample}_R2.fastq.gz` (reverse reads)

Example:
```
S1_R1.fastq.gz
S1_R2.fastq.gz
C41_S3_L004_R1_001.fastq.gz
C41_S3_L004_R2_001.fastq.gz
```

---

## **System Requirements**

### **Minimum Hardware**
- **RAM**: 64 GB (128 GB recommended)
- **Storage**: 10 TB free space
- **CPU**: 32+ cores
- **GPU**: Optional (for some deep learning applications)

### **Software Requirements**
- **OS**: Linux (Ubuntu 18.04+ or CentOS 7+)
- **Python**: 3.7-3.9
- **R**: 4.0.0+
- **Java**: 8+
- **Docker**: 20.10+

---

## **Troubleshooting**

### **Common Issues**

#### **1. Memory Errors**
```bash
# Reduce memory usage in scripts
export OMP_NUM_THREADS=4
ulimit -v 50000000  # Limit virtual memory
```

#### **2. Docker Permission Issues**
```bash
# Add user to docker group
sudo usermod -aG docker $USER
newgrp docker
```

#### **3. Reference Genome Issues**
```bash
# Verify file integrity
md5sum mm10.fa
# Re-index if corrupted
bwa index mm10.fa
samtools faidx mm10.fa
```

#### **4. R Package Installation Issues**
```r
# Install system dependencies (Ubuntu)
system("sudo apt-get install libxml2-dev libcurl4-openssl-dev libssl-dev")
# Then reinstall packages
BiocManager::install("ChIPseeker", force = TRUE)
```

### **Error Logs**

Check these log files for debugging:
- `output/{sample}/QC_check/*.log`
- `output/{sample}/STAR/*Log.final.out`
- `{sample}/outs/web_summary.html` (Cell Ranger)

---

## **Output Interpretation**

### **Key Output Files**

| Analysis Type | Key Outputs | Description |
|---------------|-------------|-------------|
| **QC** | `*_fastqc.html` | Quality reports |
| **ATAC-seq** | `*_peaks.narrowPeak` | Peak calls |
| **Single-cell** | `*.rds` | Seurat objects |
| **Hi-C** | `*_tagdir/` | Interaction data |
| **RNA-seq** | `featurecounts_output.txt` | Count matrix |

### **Visualization**

The pipeline generates various plots:
- UMAP/t-SNE plots (single-cell)
- Feature plots (gene expression)
- Quality control plots
- Peak annotation plots

---

## **Citation**

If you use this pipeline, please cite the following tools:

```bibtex
@article{seurat2021,
  title={Integrated analysis of multimodal single-cell data},
  author={Hao, Yuhan and others},
  journal={Cell},
  year={2021}
}

@article{tobias2020,
  title={TOBIAS: analysis of transcription factor binding and chromatin accessibility},
  author={Bentsen, Mette and others},
  journal={Nature Communications},
  year={2020}
}
```

---

## **Support**

### **Getting Help**

1. **Check logs**: Review error messages in output directories
2. **Documentation**: Consult individual tool documentation
3. **Issues**: Create GitHub issues with:
   - Error messages
   - System information
   - Input file examples

### **Useful Resources**

- [Seurat Documentation](https://satijalab.org/seurat/)
- [TOBIAS Documentation](https://github.com/loosolab/TOBIAS)
- [Cell Ranger Documentation](https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/what-is-cell-ranger)
- [Homer Documentation](http://homer.ucsd.edu/homer/)

---

## **Contributing**

We welcome contributions! Please:

1. Fork the repository
2. Create a feature branch
3. Submit a pull request with detailed description

### **Development Setup**

```bash
# Development environment
conda env create -f environment-dev.yml
conda activate multiomics-dev

# Run tests
bash test_pipeline.sh
```

---

## **License**

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## **Acknowledgments**

- 10x Genomics for Cell Ranger
- The Seurat team for single-cell analysis tools
- TOBIAS developers for ATAC-seq analysis
- Homer developers for Hi-C analysis tools
