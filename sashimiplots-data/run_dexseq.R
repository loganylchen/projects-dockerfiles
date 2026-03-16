#!/usr/bin/env Rscript
# ============================================
# DEXSeq Analysis Script
# Differential exon usage analysis
# ============================================

suppressPackageStartupMessages({
  library(DEXSeq)
  library(GenomicFeatures)
  library(BiocParallel)
})

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
  cat("Usage: Rscript run_dexseq.R <gtf_file> <bam_dir> <output_dir> <sample_info_file> [threads]\n")
  cat("\n")
  cat("Arguments:\n")
  cat("  gtf_file        : Path to GTF annotation file\n")
  cat("  bam_dir         : Directory containing sorted BAM files\n")
  cat("  output_dir      : Output directory for DEXSeq results\n")
  cat("  sample_info_file: Tab-separated file with columns: sample_id, condition, bam_file\n")
  cat("  threads         : Number of threads (default: 4)\n")
  quit(status = 1)
}

gtf_file <- args[1]
bam_dir <- args[2]
output_dir <- args[3]
sample_info_file <- args[4]
threads <- ifelse(length(args) >= 5, as.integer(args[5]), 4)

# Create output directory
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

cat("============================================\n")
cat("DEXSeq Differential Exon Usage Analysis\n")
cat("============================================\n")
cat("GTF file:", gtf_file, "\n")
cat("BAM directory:", bam_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("Sample info:", sample_info_file, "\n")
cat("Threads:", threads, "\n")
cat("\n")

# Set up parallel processing
BPPARAM <- MulticoreParam(workers = threads)

# Read sample information
sample_info <- read.delim(sample_info_file, stringsAsFactors = FALSE)
cat("Samples loaded:", nrow(sample_info), "\n")
print(sample_info)

# Step 1: Prepare annotation
cat("\n[Step 1] Preparing DEXSeq annotation from GTF...\n")
gff_file <- file.path(output_dir, "dexseq_annotation.gff")

# Use DEXSeq's Python script to flatten GTF
python_script <- system.file("python_scripts", "dexseq_prepare_annotation.py", 
                              package = "DEXSeq")

if (file.exists(python_script)) {
  cmd <- sprintf("python3 %s -r no %s %s", python_script, gtf_file, gff_file)
  cat("Running:", cmd, "\n")
  system(cmd)
} else {
  # Alternative: create annotation directly in R
  cat("Creating annotation using GenomicFeatures...\n")
  txdb <- makeTxDbFromGFF(gtf_file, format = "gtf")
  exons_by_gene <- exonsBy(txdb, by = "gene")
  
  # Flatten exons
  flat_exons <- disjoin(exons_by_gene)
  
  # Export as GFF
  export(flat_exons, gff_file, format = "gff3")
}

# Step 2: Count reads per exon
cat("\n[Step 2] Counting reads per exon...\n")

# Get BAM files
bam_files <- file.path(bam_dir, sample_info$bam_file)
names(bam_files) <- sample_info$sample_id

# Check BAM files exist
for (bf in bam_files) {
  if (!file.exists(bf)) {
    stop("BAM file not found: ", bf)
  }
}

# Use DEXSeq's Python script for counting
python_count_script <- system.file("python_scripts", "dexseq_count.py", 
                                    package = "DEXSeq")

count_files <- c()
for (i in seq_len(nrow(sample_info))) {
  sample_id <- sample_info$sample_id[i]
  bam_file <- bam_files[sample_id]
  count_file <- file.path(output_dir, paste0(sample_id, "_counts.txt"))
  count_files <- c(count_files, count_file)
  
  if (!file.exists(count_file)) {
    cat("Counting reads for", sample_id, "...\n")
    cmd <- sprintf("python3 %s -p yes -s no -r pos %s %s %s",
                   python_count_script, gff_file, bam_file, count_file)
    system(cmd)
  } else {
    cat("Count file exists for", sample_id, ", skipping...\n")
  }
}

# Step 3: Create DEXSeqDataSet
cat("\n[Step 3] Creating DEXSeqDataSet...\n")

# Read counts
sample_table <- data.frame(
  row.names = sample_info$sample_id,
  condition = factor(sample_info$condition)
)

dxd <- DEXSeqDataSetFromHTSeq(
  countfiles = count_files,
  sampleData = sample_table,
  design = ~ sample + exon + condition:exon,
  flattenedfile = gff_file
)

cat("DEXSeqDataSet created:\n")
cat("  Genes:", length(unique(geneIDs(dxd))), "\n")
cat("  Exons:", nrow(dxd), "\n")
cat("  Samples:", ncol(dxd), "\n")

# Step 4: Run DEXSeq analysis
cat("\n[Step 4] Running DEXSeq analysis...\n")

# Estimate size factors
cat("  Estimating size factors...\n")
dxd <- estimateSizeFactors(dxd)

# Estimate dispersions
cat("  Estimating dispersions...\n")
dxd <- estimateDispersions(dxd, BPPARAM = BPPARAM)

# Test for differential exon usage
cat("  Testing for differential exon usage...\n")
dxd <- testForDEU(dxd, BPPARAM = BPPARAM)

# Estimate exon fold changes
cat("  Estimating exon fold changes...\n")
dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "condition", BPPARAM = BPPARAM)

# Step 5: Extract results
cat("\n[Step 5] Extracting results...\n")

results <- DEXSeqResults(dxd)

# Convert to data frame
results_df <- as.data.frame(results)
results_df$geneID <- geneIDs(dxd)
results_df$exonID <- featureIDs(dxd)

# Add genomic coordinates
results_df$chr <- as.character(seqnames(results))
results_df$start <- start(results)
results_df$end <- end(results)
results_df$strand <- as.character(strand(results))

# Save results
output_file <- file.path(output_dir, "dexseq_results.csv")
write.csv(results_df, output_file, row.names = FALSE)
cat("Results saved to:", output_file, "\n")

# Save significant results
sig_results <- results_df[!is.na(results_df$padj) & results_df$padj < 0.05, ]
sig_output_file <- file.path(output_dir, "dexseq_significant.csv")
write.csv(sig_results, sig_output_file, row.names = FALSE)
cat("Significant results (FDR < 0.05):", nrow(sig_results), "\n")
cat("Saved to:", sig_output_file, "\n")

# Save DEXSeqDataSet object
rds_file <- file.path(output_dir, "dexseq_dxd.rds")
saveRDS(dxd, rds_file)
cat("DEXSeqDataSet saved to:", rds_file, "\n")

# Generate HTML report (optional)
cat("\n[Step 6] Generating HTML report...\n")
report_dir <- file.path(output_dir, "report")
tryCatch({
  DEXSeqHTML(dxd, FDR = 0.05, path = report_dir, 
             file = "DEXSeq_report.html",
             BPPARAM = BPPARAM)
  cat("HTML report generated in:", report_dir, "\n")
}, error = function(e) {
  cat("Warning: Could not generate HTML report:", conditionMessage(e), "\n")
})

cat("\n============================================\n")
cat("DEXSeq analysis complete!\n")
cat("============================================\n")
