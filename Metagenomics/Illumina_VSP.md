# Illumina Metagenomics Pipeline for Viral Surveillance Project (VSP) Data

This document describes a bioinformatics pipeline for analyzing Illumina sequencing data as part of a Viral Surveillance Project (VSP). The pipeline performs quality control, host genome filtering, and metagenomic classification using tools like **fastp**, **Bowtie2**, and **Kaiju**.

---

## Pipeline Overview

1. **Quality Control**: Trimming and filtering of raw sequencing reads using `fastp`.
2. **Host Genome Removal**: Alignment of reads to the host genome (e.g., human) using `Bowtie2` to filter out host sequences.
3. **Metagenomic Classification**: Classification of non-host reads using `Kaiju` to identify viral and microbial sequences.

---

## Prerequisites

### Tools Required

- **fastp**: For trimming and quality control of raw sequencing reads.
- **Bowtie2**: For aligning reads to a host genome and filtering out host sequences.
- **Samtools**: For handling SAM/BAM files produced by `Bowtie2`.
- **Kaiju**: For metagenomic classification of non-host reads.
- **SPAdes**: For genome assembly (optional).

### Input Data

- Raw Illumina sequencing reads in FASTQ format.
- Host reference genome (e.g., human genome) in FASTA format.

---

## Pipeline Workflow

### 1. Directory Setup

The pipeline begins by creating necessary directories to store intermediate and output files.

```bash
# Define directories
WORKING_DIR=$(pwd)
SCRIPTS="$HOME/Scripts"
CLEAN_READS="$WORKING_DIR/cleanReads"
FASTP_OUT_DIR="$WORKING_DIR/fastp_out"
KAIJU_DIR="$WORKING_DIR/kaiju_dir"
DBs="/path/to/kaijudb"  # Replace with actual path to Kaiju databases

# Create directories if they don't exist
mkdir -p "$FASTP_OUT_DIR"
mkdir -p "$CLEAN_READS"
mkdir -p "$KAIJU_DIR"
```
### 2. Quality Control with ```fastp```

Trim adapters and filter low-quality bases using fastp for high-quality data processing.

```bash
# Run fastp for quality control

$SCRIPTS/run_fastp.sh" -i "$reads_dir" -o "$FASTP_OUT_DIR" -c "$threads"
```
- ```-i``` specifies the input reads directory.
- ```-o``` specifies the output directory for cleaned reads.
- ```-c``` specifies the number of threads to use.

After running fastp, move the cleaned reads to a designated directory:

```bash
# Move and compress cleaned reads
mv "$FASTP_OUT_DIR"/*/*.fastq "$CLEAN_READS"
gzip "$CLEAN_READS"/*.fastq
chmod +rwx "$CLEAN_READS"/*.fastq
```
### 3. Host Genome Removal with Bowtie2

Remove reads that map to the host genome using Bowtie2 to ensure that only non-host reads remain for metagenomic analysis.

```bash
# Build Bowtie2 index for host genome
bowtie2-build "$reference_genome" "$reference_genome_filename" --threads "$threads"

# Create output directory for non-host reads
mkdir -p "nonHost"

# Map reads to host genome and extract unmapped reads
for fwd_file in "$CLEAN_READS"/*.fastp_1.fastq.gz; do
  base=$(basename "$fwd_file" .fastp_1.fastq.gz)
  rev_file="$CLEAN_READS/${base}.fastp_2.fastq.gz"
  
  # Output files
  sam_output="${base}/${base}.sam"
  unmapped_output="nonHost/${base}_reads_unmapped.fastq"
  
  # Run Bowtie2
  bowtie2 -1 "$fwd_file" -2 "$rev_file" -S "$sam_output" --un-conc "$unmapped_output" --threads "$threads" -x "$reference_genome_filename"
  
  echo "Mapping completed for $base"
  
  # Print mapping statistics
  samtools flagstat -@ "$threads" "$base/${base}.sam" > "$base/${base}.flagstat"
done
```
