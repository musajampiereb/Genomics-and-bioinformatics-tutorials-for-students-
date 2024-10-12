# Taxonomic classification on Illumina Metagenomics Pipeline for Viral Surveillance Pannel (VSP) Data

The Viral Surveillance Panel v2 Kit enables whole-genome sequencing (WGS) of a broad range of viruses, including SARS-CoV-2, influenza, arboviruses, and hepatitis, using a hybrid-capture method. This approach allows for viral genome sequencing with less read depth compared to shotgun metagenomics, making it ideal for viral evolution studies and broad surveillance. Unlike amplicon sequencing, hybrid-capture better identifies mutations, making it useful for monitoring rapidly evolving viruses during outbreaks. The kit integrates library preparation, target enrichment, sequencing, and data analysis, and is capable of sequencing over 200 viral pathogens. It is particularly useful for outbreak analysis, such as for Marburg virus, when sequencing is targeted at a specific pathogen.

To run the associated bioinformatics workflow, Linux or Windows Subsystem for Linux (WSL) is required. The system should have at least a Core i7 processor, 16GB of RAM, >2GHz speed, and 500GB of storage

This document describes a bioinformatics pipeline for analyzing Illumina sequencing data as part of a Viral Surveillance Pannel (VSP). The pipeline performs quality control, host genome filtering, and metagenomic classification using tools like **fastp**, **Bowtie2**, and **Kaiju**.

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
SCRIPTS="$WORKING_DIR/Scripts"
CLEAN_READS="$WORKING_DIR/cleanReads"
FASTP_OUT_DIR="$WORKING_DIR/fastp_out"
KAIJU_DIR="$WORKING_DIR/kaiju_dir"
READS_DIR="$WORKING_DIR/reads_dir"
DBs="$WORKING_DIR/databases/kaiju_viral_db"  # Replace with actual path to Kaiju databases

# Create directories if they don't exist
mkdir -p "$FASTP_OUT_DIR"
mkdir -p "$CLEAN_READS"
mkdir -p "$KAIJU_DIR"
mkdir -p "$READS_DIR"
mkdir -p "$DBs"
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
### 4. DeNovo Metagenomic Assembly

Use either ```Megahit``` , ```SPAdes``` or ```Velvet1``` (depending on the available memory) to assemble the non-host reads into contigs.

```bash
# Assemble non-host reads with Megahit
megahit -1 "nonHost/${base}_reads_unmapped.1.fastq" \
        -2 "nonHost/${base}_reads_unmapped.2.fastq" \
        -t "$threads" -o "$base/megahit_out"

# Alternatively, assemble with SPAdes
spades.py -o "$base" --meta \
          -1 "nonHost/${base}_reads_unmapped.1.fastq" \
          -2 "nonHost/${base}_reads_unmapped.2.fastq" \
          --only-assembler -t "$threads"

## When working with a limited memory, please try velvet

# Step 1: Prepare dataset with velveth
velveth "$base/velvet_out" 120 -shortPaired -fastq -separate \
    "nonHost/${base}_reads_unmapped.1.fastq" \
    "nonHost/${base}_reads_unmapped.2.fastq"

# Step 2: Perform assembly with velvetg
velvetg "$base/velvet_out" -exp_cov auto -cov_cutoff auto -ins_length 300 -min_contig_lgth 200
```
### 5. Taxonomic Classification with Kaiju

Run ```Kaiju``` to classify the assembled contigs into taxonomic categories based on the reference database.

```bash
# Run Kaiju for taxonomic classification
kaiju -t "$DBs/nodes.dmp" \
      -f "$DBs/kaiju_db_viruses.fmi" \
      -i "${base}/velvet_out/contigs.fa" -z "$threads" \
      -o "$base/${base}_kaiju.out"
```
Add taxon names to the Kaiju output:

```bash
# Add taxon names to Kaiju output
kaiju-addTaxonNames -t "$DBs/nodes.dmp" \
                    -n "$DBs/names.dmp" \
                    -i "$base/${base}_kaiju.out" \
                    -o "$base/${base}_kaiju-names.out" -p
```

Filter classified sequences:

```bash
# Extract classified sequences from Kaiju output
grep '^C' "$base/${base}_kaiju-names.out" | sed 's/;/\t/g' > "$base/${base}_classified-kaiju.tmp"
awk -v base="$base" '{OFS="\t"; print base, $0}' "$base/${base}_classified-kaiju.tmp" > "$base/${base}_classified-kaiju.out"
rm "$base/${base}_classified-kaiju.tmp"
```
### 6. Compile Results

Finally, compile the results from all samples into a single file for further analysis.

```bash
# Compile results across all samples
cat */*_classified-kaiju.out | awk '$2 == "C"' > allSamples_kaiju_results.txt

echo "Metagenomics analysis pipeline completed!"
```
