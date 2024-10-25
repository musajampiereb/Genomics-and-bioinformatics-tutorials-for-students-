# Read-level Taxonomic classification on Illumina Viral Surveillance Pannel (VSP) Data

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
chmod +rwx "$CLEAN_READS"/*.fastq.gz
```
### 3. Host Genome Removal (dehost) with ```Bowtie2```

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
### 4. Taxonomic Classification with Kaiju

Since we are performing read-level classification, there is no need to assemble the reads into contigs. Instead, classify the non-host reads directly using ```Kaiju```.

```bash
#!/bin/bash
# Define paths
DB_PATH="$DBs/kaiju_db_viruses.fmi"  # Path to Kaiju virus database
NODES_PATH="$DBs/nodes.dmp"          # Path to Kaiju taxonomic nodes file
NAMES_PATH="$DBs/names.dmp"          # Path to Kaiju taxonomic names file
THREADS=12

# Perform classification on each de-hosted sample
for fwd_file in "$WORKING_DIR/reference_genomes/human/nonHost"/*_reads_unmapped.1.fastq; do
    base=$(basename "$fwd_file" _reads_unmapped.1.fastq)
    rev_file="$WORKING_DIR/reference_genomes/human/nonHost/${base}_reads_unmapped.2.fastq"

    # Run Kaiju for taxonomic classification of de-hosted reads
    kaiju -z "$THREADS" -t "$NODES_PATH" -f "$DB_PATH" -i "$fwd_file" -j "$rev_file" -o "$KAIJU_DIR/${base}_kaiju.out"

    # Add taxon names to the Kaiju output
    kaiju-addTaxonNames -t "$NODES_PATH" -n "$NAMES_PATH" -i "$KAIJU_DIR/${base}_kaiju.out" -o "$KAIJU_DIR/${base}_kaiju-names.out"
    
    echo "Kaiju classification completed for $base"
done
```

### 5. Filter and Extract Virus Classifications

After classification, extract viral sequences from the Kaiju output. This can be done by filtering for the viral taxonomic groups.

```bash
# Filter Kaiju output for viral classifications
for file in "$KAIJU_DIR"/*_kaiju-names.out; do
    base=$(basename "$file" _kaiju-names.out)

    # Display the contents of the input file (first 20 lines)
    echo "Processing file: $file"
    echo "Contents of $file (first 20 lines):"
    head -n 20 "$file"

    # Try extracting lines classified as viruses using different patterns
    echo "Searching for virus classifications..."
    grep -i -E 'virus|viruses' "$file" > "$KAIJU_DIR/${base}_viruses.out"

    # Output the number of lines found
    found_lines=$(wc -l < "$KAIJU_DIR/${base}_viruses.out")
    echo "Lines found for viruses: $found_lines"

    # Check if output file is empty
    if [ ! -s "$KAIJU_DIR/${base}_viruses.out" ]; then
        echo "No virus classifications found for $base."
    else
        echo "Virus classification extracted for $base"
    fi

    # Optional: Show the first few lines of the output file
    echo "First few lines of ${base}_viruses.out:"
    head -n 5 "$KAIJU_DIR/${base}_viruses.out"
    echo "-----------------------------"
done

```
### 6. Relative Abundance Calculation and Stacked Bar Plot

For each sample, calculate the relative abundance of each viral taxon based on the number of classified reads.

```
# Calculate relative abundance of viruses
for file in "$KAIJU_DIR"/*_viruses.out; do
    base=$(basename "$file" _viruses.out)

    # Count the total number of reads classified as viruses
    total_virus_reads=$(wc -l < "$file")

    # Count the reads for each viral taxon and calculate relative abundance
    awk -F'\t' '{print $3}' "$file" | sort | uniq -c | awk -v total="$total_virus_reads" '{OFS="\t"; print $2, $1, ($1/total)*100}' > "$KAIJU_DIR/${base}_virus_abundance.txt"

    echo "Relative abundance calculated for $base"
done
```

### 7. Plotting the Stacked Bar Chart (R Script)

Use the ```ggplot2``` library in R to generate a stacked bar chart showing the relative abundance of viral taxa for each sample.

```
# R Script for Plotting Relative Abundance

# Load required library
library(ggplot2)

# Define the path to the abundance data file
abundance_file_path <- "/Volumes/Biospace/MVD_Rwanda_VSP/MVD_assemblies/kaiju_dir/Sample016_S16_virus_abundance.txt"

# Check if the file exists
if (!file.exists(abundance_file_path)) {
    stop("Error: The abundance file does not exist at the specified path.")
}

# Read the data
abundance_data <- read.table(abundance_file_path, header = TRUE, sep = "\t")

# Check if the data is empty
if (nrow(abundance_data) == 0) {
    stop("Error: The abundance data file is empty.")
}

# Check for required columns
required_columns <- c("Sample", "Abundance", "Taxon")
if (!all(required_columns %in% colnames(abundance_data))) {
    stop("Error: The abundance data file must contain the following columns: Sample, Abundance, Taxon.")
}

# Create a stacked bar chart
ggplot(abundance_data, aes(x = Sample, y = Abundance, fill = Taxon)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    labs(title = "Viral Taxa Distribution",
         x = "Sample",
         y = "Relative Abundance (%)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot to a file (optional)
ggsave("viral_taxa_distribution_Sample016.png", width = 10, height = 6)  # Adjust dimensions as needed

```
