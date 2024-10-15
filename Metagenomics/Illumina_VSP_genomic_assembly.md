# Reference-based genome assembly for pathogens sequenced by Illumina Viral Surveillance Pannel (VSP)

The Viral Surveillance Panel v2 Kit enables whole-genome sequencing (WGS) of a broad range of viruses, including SARS-CoV-2, influenza, arboviruses, and hepatitis, using a hybrid-capture method. This approach allows for viral genome sequencing with less read depth compared to shotgun metagenomics, making it ideal for viral evolution studies and broad surveillance. Unlike amplicon sequencing, hybrid-capture better identifies mutations, making it useful for monitoring rapidly evolving viruses during outbreaks. The kit integrates library preparation, target enrichment, sequencing, and data analysis, and is capable of sequencing over 200 viral pathogens. It is particularly useful for outbreak analysis, such as for Marburg virus, when sequencing is targeted at a specific pathogen.

To run the associated bioinformatics workflow, Linux or Windows Subsystem for Linux (WSL) is required. The system should have at least a Core i7 processor, 16GB of RAM, >2GHz speed, and 500GB of storage

This document describes a bioinformatics pipeline for analyzing Illumina sequencing data as part of a Viral Surveillance Pannel (VSP). The pipeline performs quality control, host genome filtering, and metagenomic classification using tools like **fastp** and **Bowtie2**.

---

## Pipeline Overview

1. **Quality Control**: Trimming and filtering of raw sequencing reads using `fastp`.
2. **Host Genome Removal**: Alignment of reads to the host genome (e.g., human) using `Bowtie2` to filter out host sequences.
3. **Genome assembly**: Classification of non-host reads using `Kaiju` to identify viral and microbial sequences.

---

## Prerequisites

### Tools Required

- **fastp**: For trimming and quality control of raw sequencing reads.
- **Bowtie2**: For aligning reads to a host genome and filtering out host sequences.
- **Samtools**: For handling SAM/BAM files produced by `Bowtie2`.
- **minimap2**: For aligning the unmapped reads to the reference.
- **picard**: For removing duplicates.
- **LoFred**: For variant calling.
- **ivar**: For variant calling and consensus sequence generation.
        
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
SCRIPTS="$WORKING_DIR/Scripts"
CLEAN_READS="$WORKING_DIR/cleanReads"
FASTP_OUT_DIR="$WORKING_DIR/fastp_out"
KAIJU_DIR="$WORKING_DIR/kaiju_dir"
READS_DIR="$WORKING_DIR/reads_dir"
DBs="$WORKING_DIR/databases/kaiju_viral_db"  # Replace with actual path to Kaiju databases

# Create directories if they don't exist
mkdir -p "$FASTP_OUT_DIR"
mkdir -p "$CLEAN_READS"
mkdir -p "$READS_DIR"
mkdir -p "$SCRIPTS"
```
### 2. Quality Control with ```fastp```

Trim adapters and filter low-quality bases using ```fastp``` for high-quality data processing.

```bash
# Run fastp for quality control
threads=4
"$SCRIPTS/run_fastp.sh" -i "$READS_DIR" -o "$FASTP_OUT_DIR" -c "$threads"
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
### 3. Host Genome Removal with ```Bowtie2```

Remove reads that map to the host genome using ```Bowtie2``` to ensure that only non-host reads remain for metagenomic analysis.

```bash
#!/bin/bash

# Step 1: Create directories for reference genomes and download the human reference genome
mkdir -p $WORKING_DIR/reference_genomes/human
cd $WORKING_DIR/reference_genomes/human

# Download the human reference genome
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

# Unzip the genome
gunzip hg38.fa.gz

# Set reference genome path variables
reference_genome="$WORKING_DIR/reference_genomes/human/hg38.fa.gz"
index_prefix="$WORKING_DIR/reference_genomes/human/index"

# Step 2: Build Bowtie2 index for host genome
bowtie2-build "$reference_genome" "$index_prefix" --threads "$threads"

# Step 3: Create output directory for non-host reads
mkdir -p "$WORKING_DIR/reference_genomes/human/nonHost"

# Step 4: Map reads to host genome and extract unmapped reads
for fwd_file in "$CLEAN_READS"/*.fastp_1.fastq.gz; do
  base=$(basename "$fwd_file" .fastp_1.fastq.gz)
  rev_file="$CLEAN_READS/${base}.fastp_2.fastq.gz"

  # Create output directory for each sample
  mkdir -p "$WORKING_DIR/reference_genomes/human/${base}"

  # Output files
  sam_output="${base}/${base}.sam"
  unmapped_output="nonHost/${base}_reads_unmapped.fastq"

  # Step 5: Run Bowtie2
  bowtie2 -1 "$fwd_file" -2 "$rev_file" -S "$WORKING_DIR/reference_genomes/human/$sam_output" --un-conc "$WORKING_DIR/reference_genomes/human/$unmapped_output" --threads 12 -x "$index_prefix"

  echo "Mapping completed for $base"

  # Step 6: Print mapping statistics
  samtools flagstat -@ "$threads" "$WORKING_DIR/reference_genomes/human/$sam_output" > "$WORKING_DIR/reference_genomes/human/${base}/${base}.flagstat"
done
```

### 4. Mapping non-human reads to the reference sequence

In this step the de-hosted reads are mapped to the reference sequence. Based on the sequence run or target virus the reference sequence can be downloaded from the NCBI-Virus database. The information on which reference (accession number) to be used for the selected virus can be searched in literature normally has a refseq tag in the NCBI search results.

### Mapping using minimap2 
Here is the syntax for mapping non-human reads to the reference sequence using minimap2:

```
# Define the reference genome:
reference="$WORKING_DIR/reference_genomes/human/NC_001608.3"
r1="nonHost/${base}_reads_unmapped.1.fastq"
r2="nonHost/${base}_reads_unmapped.2.fastq"

# Define aligned SAM file output:
outfile="$WORKING_DIR/reference_genomes/human/nonHost/aligned.sam"

# Run Minimap2 command:
minimap2 -ax sr ${reference} ${r1} ${r2} > ${outfile}

```

### 5. Processing of the mapped sequences

### Sort SAM file
```
infile="$WORKING_DIR/reference_genomes/human/nonHost/aligned.sam"
outfile="$WORKING_DIR/reference_genomes/human/nonHost/aligned.sorted.bam"
samtools sort ${infile} > ${outfile}
```
### Discard un-mapped reads
```
infile="$WORKING_DIR/reference_genomes/human/nonHost/aligned.sorted.bam"
outfile="$WORKING_DIR/reference_genomes/human/nonHost/aligned.sorted.mapped.bam"
samtools view -F 0x04 -b ${infile} > ${outfile}
```
### Index BAM file
```
samtools index ${outfile}

#In case you install Picard through conda:

conda activate picard

#Tag duplicate reads in BAM file:
infile="aligned.sorted.mapped.bam"
outfile="aligned.sorted.mapped.markduplicates.bam"
outmetrics="aligned.sorted.mapped.markduplicates.metrics.txt"
picard MarkDuplicates \
 -Xmx8g \
 -I ${infile} \
 -O ${outfile} \
 -M ${outmetrics}
# Index BAM file:
samtools index ${outfile}
```
### 6. Variant calling

Variant calling is the process of identifying and cataloging the differences between the virus of interest sequencing reads and a reference genome. Variant calling enables us to know the amount of changes occurred on the genome of interest, we get the SNPs. MNPs, idels e.t.c. This step is crucial for knowing/detection of new variants of the virus in question.

#### 6.1 variant calling using ivar

In case you install iVar through conda environment (ivar_env):
```
conda activate ivar_env
```
Make a pileup and pipe to iVar to call variants:
```
infile="aligned.sorted.mapped.markduplicates.bam"
prefix="out_variants"
samtools mpileup --reference ${reference} ${infile} | ivar variants -r ${reference} -p ${prefix}
```
#### 6.2 variant calling using Lofreq

In case you install LoFreq through conda environment:
```
conda activate lofreq_env
```
Call variants
infile="aligned.sorted.mapped.bam"
outfile="variants.vcf"
lofreq call -f ${reference} -o ${outfile} ${infile}

### 7. Consensus Calling

At this step we now create the consensus fasta file which we will use for downstream analysis, here we call consensus for position with support from at least on read at that position, this is because we use metagenomic sequencing aproach the depth might be very low. if tilling approach was used for sequencing depth for calling consensus can be set to 5 or 10 reads per position.

```
# In case you install iVar through conda:
conda activate ivar

# Generate consensus FASTA:
# Optionally, you can set different parameters to define minimum thresholds for the consensus

# (see ivar consensus help).
infile="aligned.sorted.mapped.bam"
outfile="consensus_sequence.fa"
samtools mpileup -A -Q 0 ${infile} | ivar consensus -p ${outfile} -q 10 -t 0 -m 1

```
