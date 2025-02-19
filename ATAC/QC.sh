#######################
#$1="sample name: example "sample_S1"
#$2="work directory"
######################
#!/usr/bin/env bash
set -o errexit
set -o nounset
mkdir -p output
mkdir -p output/${1}
#fastqc
mkdir -p output/${1}/QC_check
fastqc -t 30 ${1}_R1.fastq.gz -o output/${1}/QC_check
fastqc -t 30 ${1}_R2.fastq.gz -o output/${1}/QC_check
mkdir -p output/$1/trimedreads_fastp
#sudo chmod 777 output/$1/trimedreads_fastp
fastp -i ${1}_R1.fastq.gz -I ${1}_R2.fastq.gz -o output/$1/trimedreads_fastp/${1}_1.fq.gz -O output/$1/trimedreads_fastp/${1}_2.fq.gz -j output/$1/trimedreads_fastp/${1}.json -h output/$1/trimedreads_fastp/${1}.html -w 16
#Fastqc_fastp
mkdir -p output/${1}/Fastp_fastqc
fastqc -t 30 output/$1/trimedreads_fastp/${1}_1.fq.gz -o output/${1}/Fastp_fastqc
fastqc -t 30 output/$1/trimedreads_fastp/${1}_2.fq.gz -o output/${1}/Fastp_fastqc

