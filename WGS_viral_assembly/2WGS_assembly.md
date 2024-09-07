# Whole Genome Sequence Assembly Workflow for Viral Genome (MPOX)

**NOTE1:** This assignment is based on published raw data from the MPOX outbreak in Kamituga, South Kivu.  
Access the publication here: [https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2024.29.11.2400106](https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2024.29.11.2400106)  

**NOTE2:** This workflow requires "conda" installed and running on your computer.  

This workflow guides you through downloading raw sequencing data, performing quality control, trimming, alignment, variant calling, and finally creating a consensus sequence for the MPOX virus.

## Step 0: Create and Activate Environment, Set Up Working Directories, and Environment Variables

1. Create conda environment

    ```bash
    conda create -n wgs_env python=3.8
    ```

2. Activate the created conda environment

    ```bash
    conda activate wgs_env
    ```

3. Create a working directory and set the `WGS_HOME` environment variable

    ```bash
    mkdir -p ~/workspace/bioinformatics/
    export WGS_HOME=~/workspace/bioinformatics
    ```

4. Ensure that the working directory is set correctly

    ```bash
    echo "Working Directory: $WGS_HOME"
    ```

5. Set up additional environment variables for different directories in the workflow

    ```bash
    export WGS_DATA_DIR=$WGS_HOME/data
    export WGS_DATA_TRIM_DIR=$WGS_DATA_DIR/trimmed
    export WGS_REFS_DIR=$WGS_HOME/refs
    export WGS_REF_FASTA=$WGS_REFS_DIR/mpox_ref.fa
    export WGS_ALIGN_DIR=$WGS_HOME/alignments/minimap2
    export WGS_DEPTH_DIR=$WGS_HOME/depth
    export WGS_VCF_DIR=$WGS_HOME/vcf
    ```

6. Ensure that all environment variables are correctly defined

    ```bash
    env | grep WGS
    ```

## Step 1: Download Raw Data and Reference Genome

1. Create the necessary directory for raw data and navigate to it

    ```bash
    mkdir -p $WGS_DATA_DIR
    cd $WGS_DATA_DIR
    ```

2. Download the raw FASTQ data files

    ```bash
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR126/001/ERR12670101/ERR12670101.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR126/004/ERR12670104/ERR12670104.fastq.gz
    wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR126/005/ERR12670105/ERR12670105.fastq.gz
    ```

3. Unzip the downloaded FASTQ files for inspection

    ```bash
    gunzip *.fastq.gz
    ```

4. Rename the FASTQ file for easier reference in subsequent steps

    ```bash
    mv ERR12670104.fastq all_reads.fastq
    ```

5. Download the reference genome for the MPOX virus

    ```bash
    cd $WGS_REFS_DIR
    wget https://www.ebi.ac.uk/ena/browser/api/fasta/JX878425.1?download=true -O mpox_ref.fa
    ```

## Step 2: Quality Control (QC) Check Using NanoPlot

1. Ensure that NanoPlot is installed. If not, install it via conda

    ```bash
    conda install bioconda::nanoplot
    ```

2. Perform a QC check on the FASTQ file using NanoPlot

    ```bash
    cd $WGS_DATA_DIR
    gzip all_reads.fastq
    NanoPlot --fastq all_reads.fastq.gz -o QC_REPORT --plots kde
    open QC_REPORT
    ```

## Step 3: Data Filtering Based on Quality Score Using fastp

1. Unzip the FASTQ file

    ```bash
    gunzip all_reads.fastq.gz
    ```

2. Install `fastp` via conda

    ```bash
    conda install bioconda::fastp
    ```

3. Create a directory for trimmed data

    ```bash
    mkdir -p $WGS_DATA_TRIM_DIR
    ```

4. Filter the FASTQ data using `fastp`

    ```bash
    fastp -w 48 -i all_reads.fastq -l 100 -q 09 -o $WGS_DATA_TRIM_DIR/all_readsQC.fastq
    ```

## Step 4: Trim Adapters Using Cutadapt

1. Navigate to the trimmed data directory

    ```bash
    cd $WGS_DATA_TRIM_DIR
    ```

2. Install `cutadapt` via conda

    ```bash
    conda install bioconda::cutadapt
    ```

3. Trim adapters from the FASTQ data using `cutadapt`

    ```bash
    cutadapt -u 30 -o all_reads_QC1.fastq all_readsQC.fastq
    cutadapt -u -30 -o all_reads_QC2.fastq all_reads_QC1.fastq
    ```

## Step 5: Align Reads to the Reference Genome Using Minimap2

1. Create a directory for alignments and navigate to it

    ```bash
    mkdir -p $WGS_ALIGN_DIR
    cd $WGS_ALIGN_DIR
    ```

2. Install `minimap2` via conda

    ```bash
    conda install bioconda::minimap2
    ```

3. Align reads to the reference genome using `minimap2`

    ```bash
    minimap2 -Y -t 12 -x map-ont -a $WGS_REF_FASTA $WGS_DATA_TRIM_DIR/all_reads_QC2.fastq > all_reads.sam
    ```

## Step 6: Sort the SAM File Using Samtools

1. Install `samtools` via conda

    ```bash
    conda install bioconda::samtools
    ```

2. Sort the SAM file

    ```bash
    samtools sort -O SAM -o all_reads_sorted.sam all_reads.sam
    ```

## Step 7: Convert SAM to BAM and Sort

1. Convert the SAM file to BAM format and sort

    ```bash
    samtools view -bS all_reads_sorted.sam | samtools sort -@ 16 -o all_reads.bam
    ```

## Step 8: Index the BAM File

1. Index the BAM file

    ```bash
    samtools index all_reads.bam
    ```

## Step 9: Create a Depth Profile

1. Create a directory for depth profiles and navigate to it

    ```bash
    mkdir -p $WGS_DEPTH_DIR
    cd $WGS_DEPTH_DIR
    ```

2. Generate a depth profile from the BAM file

    ```bash
    samtools mpileup -a -A -Q 0 -d 0 -f $WGS_REF_FASTA $WGS_ALIGN_DIR/all_reads.bam | awk '{print $2, $3, $4}' > all_reads.depth
    ```

## Step 10: Variant Calling Using BCFtools

1. Install `bcftools` via conda

    ```bash
    conda install bioconda::bcftools
    ```

2. Create a directory for VCF files and navigate to it

    ```bash
    mkdir -p $WGS_VCF_DIR
    cd $WGS_VCF_DIR
    ```

3. Call variants using `bcftools`

    ```bash
    bcftools mpileup -f $WGS_REF_FASTA $WGS_ALIGN_DIR/all_reads.bam | bcftools call -mv -Oz -o all_reads.vcf.gz
    ```

## Step 11: Filter the VCF File for High-Quality Variants

1. Filter the VCF file to retain high-quality variants

    ```bash
    bcftools filter -i 'DP>20 && AF>0.1' all_reads.vcf.gz -Oz -o all_reads_filtered.vcf.gz
    ```

## Step 12: Create a Consensus Sequence

#!/bin/bash
### Consensus Sequence Generation Script

### Define file paths

```bash
depth_file="$WGS_DEPTH_DIR/all_reads.depth"
input_fasta="all_reads_consensus_temp.fasta"
output_fasta="all_reads_consensus.fasta"
filtered_consensus="all_reads_consensus_filtered.fasta"
```
### Minimum depth threshold
min_depth=20

### Filter the depth profile and prepare regions to include

```bash
awk -v min_depth=$min_depth '
BEGIN {
    while ((getline < depth_file) > 0) {
        if ($3 >= min_depth) {
            depth[$2] = $3
        }
    }
}
{
    if (/^>/) {
        print $0
    } else {
        seq = $0
        for (i = 1; i <= length(seq); i++) {
            if (!(i in depth)) {
                seq = substr(seq, 1, i-1) "N" substr(seq, i+1)
            }
        }
        print seq
    }
}' "$input_fasta" > "$output_fasta"
```

### Install ivar if not already installed
```bash
conda install -c bioconda ivar
```
### Generate consensus sequence using ivar
```bash
ivar consensus -t 0 -i $WGS_ALIGN_DIR/all_reads.bam -f $WGS_REF_FASTA -b $WGS_VCF_DIR/all_reads_filtered.vcf.gz -o $filtered_consensus

echo "Consensus sequence generated and saved to $filtered_consensus"
```
