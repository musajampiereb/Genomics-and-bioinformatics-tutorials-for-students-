```bash
#!/bin/bash
#######################################################################
## Whole Genome Sequence Assembly Workflow for Viral Genome (MPOX)  ##
#######################################################################

# NOTE1: This assignment is based on published raw data from the MPOX outbreak in Kamituga, South Kivu.
# Access the publication here: https://www.eurosurveillance.org/content/10.2807/1560-7917.ES.2024.29.11.2400106 

# NOTE2: This workflow requires "conda" to be installed and running on your computer 

# This workflow guides you through downloading raw sequencing data, performing quality control, trimming, alignment, 
# variant calling, and finally creating a consensus sequence for the MPOX virus.

# Step 0: Create and activate Environment, set up Working Directories, and Environment Variables
###############################################################
# Create conda environment
conda create -n wgs_env python=3.8

# Activate the created conda environment 
conda activate wgs_env

# Create a working directory and set the ‘WGS_HOME’ environment variable.
mkdir -p ~/workspace/bioinformatics/
export WGS_HOME=~/workspace/bioinformatics

# Ensure that the working directory is set correctly.
echo "Working Directory: $WGS_HOME"

# Set up additional environment variables for different directories in the workflow.
export WGS_DATA_DIR=$WGS_HOME/data
export WGS_DATA_TRIM_DIR=$WGS_DATA_DIR/trimmed
export WGS_REFS_DIR=$WGS_HOME/refs
export WGS_REF_FASTA=$WGS_REFS_DIR/mpox_ref.fa
export WGS_ALIGN_DIR=$WGS_HOME/alignments/minimap2
export WGS_DEPTH_DIR=$WGS_HOME/depth
export WGS_VCF_DIR=$WGS_HOME/vcf

# Ensure that all environment variables are correctly defined.
env | grep WGS

# Step 1: Download Raw Data and Reference Genome
################################################

# Create the necessary directory for raw data and navigate to it.
mkdir -p $WGS_DATA_DIR
cd $WGS_DATA_DIR

# Download the raw FASTQ data files.
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR126/001/ERR12670101/ERR12670101.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR126/004/ERR12670104/ERR12670104.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR126/005/ERR12670105/ERR12670105.fastq.gz

# Unzip the downloaded FASTQ files for inspection.
gunzip *.fastq.gz

# Rename the FASTQ file for easier reference in subsequent steps.
mv ERR12670104.fastq all_reads.fastq

# Download the reference genome for the MPOX virus.
cd $WGS_REFS_DIR
wget https://www.ebi.ac.uk/ena/browser/api/fasta/JX878425.1?download=true -O mpox_ref.fa

# Step 2: Quality Control (QC) Check Using NanoPlot
###################################################

# Go back to the folder that stores your data.
cd $WGS_DATA_DIR

# Ensure that NanoPlot is installed. If not, install it via conda.
conda install bioconda::nanoplot

# Perform a QC check on the FASTQ file using NanoPlot.
gzip all_reads.fastq
NanoPlot --fastq all_reads.fastq.gz -o QC_REPORT --plots kde

# Step 3: Data Filtering Based on Quality Score Using fastp
###########################################################

# Ensure the FASTQ file is unzipped.
gunzip all_reads.fastq.gz

# Install FASTP tool 
conda install bioconda::fastp

# Create a directory to store quality control reports and navigate to it.
mkdir -p $WGS_DATA_TRIM_DIR

# Filter the reads based on read quality with a quality score of 9.
fastp -w 48 -i all_reads.fastq -l 100 -q 09 -o $WGS_DATA_TRIM_DIR/all_readsQC.fastq

# Step 4: Trim Adapters Using Cutadapt
#######################################

# Change the directory to where the trimmed data are stored.
cd $WGS_DATA_TRIM_DIR

# Ensure cutadapt is installed.
conda install bioconda::cutadapt

# Trim adapters from both ends of the reads.
cutadapt -u 30 -o all_reads_QC1.fastq all_readsQC.fastq
cutadapt -u -30 -o all_reads_QC2.fastq all_reads_QC1.fastq

# Step 5: Align Reads to the Reference Genome Using Minimap2
############################################################

# Create the alignment directory and navigate to it.
mkdir -p $WGS_ALIGN_DIR
cd $WGS_ALIGN_DIR

# Align the reads to the reference genome using minimap2 and generate a SAM file.
conda install bioconda::minimap2
minimap2 -Y -t 12 -x map-ont -a $WGS_REF_FASTA $WGS_DATA_TRIM_DIR/all_reads_QC2.fastq > all_reads.sam

# Step 6: Sort the SAM File Using Samtools
###########################################

# Sort the SAM file using samtools.
conda install bioconda::samtools
samtools sort -O SAM -o all_reads_sorted.sam all_reads.sam

# Step 7: Convert SAM to BAM and Sort
#####################################

# Convert the sorted SAM file to BAM format and sort the BAM file.
samtools view -bS all_reads_sorted.sam | samtools sort -@ 16 -o all_reads.bam

# Step 8: Index the BAM File
#############################

# Index the BAM file to prepare it for further analysis.
samtools index all_reads.bam

# Step 9: Create a Depth Profile
################################

# Create a directory to store depth information.
mkdir -p $WGS_DEPTH_DIR
cd $WGS_DEPTH_DIR

# Generate a depth profile of the aligned reads.
samtools mpileup -a -A -Q 0 -d 0 -f $WGS_REF_FASTA $WGS_ALIGN_DIR/all_reads.bam | awk '{print $2, $3, $4}' > all_reads.depth

# Step 10: Variant Calling Using BCFtools
#########################################

# Install bcftools with conda.
conda install bioconda::bcftools

# Create a directory to store VCF files.
mkdir -p $WGS_VCF_DIR
cd $WGS_VCF_DIR

# Call variants using bcftools.
bcftools mpileup -f $WGS_REF_FASTA $WGS_ALIGN_DIR/all_reads.bam | bcftools call -mv -Oz -o all_reads.vcf.gz

# Step 11: Filter the VCF File for High-Quality Variants
#######################################################

# Filter the VCF file for high-quality variants.
bcftools filter -i 'DP>20 && AF>0.1' all_reads.vcf.gz -Oz -o all_reads_filtered.vcf.gz

# Step 12: Create a Consensus Sequence
######################################

# Index the VCF file.
bcftools index all_reads_filtered.vcf.gz

# Generate a consensus sequence using the filtered variants.
bcftools consensus -f $WGS_REF_FASTA all_reads_filtered.vcf.gz > all_reads_consensus_temp.fasta

# Incorporate depth profile information into the consensus sequence.
depth_file="$WGS_DEPTH_DIR/all_reads.depth"
input_fasta="all_reads_consensus_temp.fasta"
output_fasta="all_reads_consensus.fasta"

awk -v depth_file="$depth_file" '
BEGIN {
    # Load the depth file into an array
    while ((getline < depth_file) > 0) {
        depth[$2] = $3
    }
}
{
    if (/^>/) {
        # Print the header lines as they are
        print $0
    } else {
        seq = $0
        for (i = 1; i <= length(seq); i++) {
            # Only replace with 'N' if the depth is less than 20
            if (depth[i] < 20) {
                seq = substr(seq, 1, i-1) "N" substr(seq, i+1)
            }
        }
        print seq
    }
}' "$input_fasta" > "$output_fasta"

# Step 13: Final Report and Interpretation
##########################################

# Summarize the workflow and analyze the generated consensus sequence.


