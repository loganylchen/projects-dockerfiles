#!/bin/bash
# ============================================
# Configuration file for RNA-seq analysis
# Modify these settings according to your data
# ============================================

# ============================================
# Sample Configuration
# 2 conditions x 3 replicates = 6 samples total
# Paired-end (PE) data
# ============================================

# Condition names
CONDITION1="Control"
CONDITION2="Treatment"

# Sample names (modify as needed)
declare -A SAMPLES
# Format: SAMPLES[sample_name]="condition"
SAMPLES["Sample1"]="${CONDITION1}"
SAMPLES["Sample2"]="${CONDITION1}"
SAMPLES["Sample3"]="${CONDITION1}"
SAMPLES["Sample4"]="${CONDITION2}"
SAMPLES["Sample5"]="${CONDITION2}"
SAMPLES["Sample6"]="${CONDITION2}"

# FASTQ file paths (paired-end)
# Modify these paths to match your data location
# R1 = Read 1, R2 = Read 2
declare -A FASTQ_R1
declare -A FASTQ_R2

FASTQ_R1["Sample1"]="/data/fastq/Sample1_R1.fastq.gz"
FASTQ_R2["Sample1"]="/data/fastq/Sample1_R2.fastq.gz"

FASTQ_R1["Sample2"]="/data/fastq/Sample2_R1.fastq.gz"
FASTQ_R2["Sample2"]="/data/fastq/Sample2_R2.fastq.gz"

FASTQ_R1["Sample3"]="/data/fastq/Sample3_R1.fastq.gz"
FASTQ_R2["Sample3"]="/data/fastq/Sample3_R2.fastq.gz"

FASTQ_R1["Sample4"]="/data/fastq/Sample4_R1.fastq.gz"
FASTQ_R2["Sample4"]="/data/fastq/Sample4_R2.fastq.gz"

FASTQ_R1["Sample5"]="/data/fastq/Sample5_R1.fastq.gz"
FASTQ_R2["Sample5"]="/data/fastq/Sample5_R2.fastq.gz"

FASTQ_R1["Sample6"]="/data/fastq/Sample6_R1.fastq.gz"
FASTQ_R2["Sample6"]="/data/fastq/Sample6_R2.fastq.gz"

# ============================================
# Required Input Files (YOU MUST PROVIDE THESE)
# ============================================

# Reference genome (FASTA)
REFERENCE_FASTA="/reference/genome.fa"

# GTF annotation file
GTF_FILE="/reference/annotation.gtf"

# HISAT2 index prefix (YOU PROVIDED THIS)
HISAT2_INDEX="/reference/hisat2_index/genome"

# ============================================
# Auto-generated Files (SCRIPT WILL CREATE THESE)
# ============================================

# Transcriptome FASTA (will be extracted from GTF + genome)
TRANSCRIPTOME_FASTA="/reference/transcriptome.fa"

# Salmon index (will be built from transcriptome)
SALMON_INDEX="/reference/salmon_index"

# GFF3 annotation file (will be converted from GTF for MISO)
GFF3_FILE="/reference/annotation.gff3"

# ============================================
# Directory Settings
# ============================================

# Input/Output directories
DATA_DIR="/data"
OUTPUT_DIR="/output"
FASTQ_DIR="${DATA_DIR}/fastq"
REFERENCE_DIR="/reference"

# Output subdirectories
TRIMMED_DIR="${OUTPUT_DIR}/trimmed"
ALIGNED_DIR="${OUTPUT_DIR}/aligned"
SALMON_DIR="${OUTPUT_DIR}/salmon"
STRINGTIE_DIR="${OUTPUT_DIR}/stringtie"
COUNTS_DIR="${OUTPUT_DIR}/counts"
QC_DIR="${OUTPUT_DIR}/qc"
JUNCTIONS_DIR="${OUTPUT_DIR}/junctions"
LOGS_DIR="${OUTPUT_DIR}/logs"

# Splicing tool output directories
RMATS_DIR="${OUTPUT_DIR}/rmats"
SUPPA2_DIR="${OUTPUT_DIR}/suppa2"
LEAFCUTTER_DIR="${OUTPUT_DIR}/leafcutter"
DEXSEQ_DIR="${OUTPUT_DIR}/dexseq"
JUNCTIONSEQ_DIR="${OUTPUT_DIR}/junctionseq"
MAJIQ_DIR="${OUTPUT_DIR}/majiq"
MISO_DIR="${OUTPUT_DIR}/miso"

# ============================================
# Analysis Parameters
# ============================================

# Number of threads
THREADS=8

# HISAT2 parameters
HISAT2_OPTS="--dta --no-mixed --no-discordant"

# Read length (for rMATS)
READ_LENGTH=150

# rMATS parameters
RMATS_OPTS="--novelSS --mil 50 --mel 500"

# Leafcutter parameters
LEAFCUTTER_MIN_COVERAGE=30
LEAFCUTTER_INTRON_MAX=500000

# SUPPA2 parameters
SUPPA2_EVENT_TYPES="SE SS MX RI FL"

# DEXSeq/JunctionSeq FDR cutoff
FDR_CUTOFF=0.05

# Minimum junction read count for output
MIN_JUNCTION_READS=5

# ============================================
# Pipeline Control
# Set to 1 to run, 0 to skip
# ============================================

RUN_QC=1
RUN_TRIMMING=1
RUN_ALIGNMENT=1
RUN_SORTING=1
RUN_SALMON=1
RUN_STRINGTIE=1
RUN_COUNTS=1
RUN_JUNCTION_EXTRACT=1

# Splicing analysis tools
RUN_RMATS=1
RUN_SUPPA2=1
RUN_LEAFCUTTER=1
RUN_DEXSEQ=1
RUN_JUNCTIONSEQ=1
RUN_MAJIQ=0      # Requires license, disabled by default
RUN_MISO=1
