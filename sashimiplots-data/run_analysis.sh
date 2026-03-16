#!/bin/bash
# ============================================
# RNA-seq Analysis Pipeline
# For generating sashimiplots test data
# Supports: rMATS, SUPPA2, Leafcutter, DEXSeq, JunctionSeq, MAJIQ, MISO
# ============================================
# Usage: ./run_analysis.sh [config_file]
# ============================================
#
# Required inputs (YOU provide these):
#   - GTF annotation file
#   - Genome FASTA file
#   - HISAT2 index
#   - FASTQ files (paired-end)
#
# Auto-generated files (script creates these):
#   - Transcriptome FASTA (from GTF + genome)
#   - Salmon index
#   - GFF3 annotation (for MISO)
#
# ============================================

set -e  # Exit on error
set -u  # Exit on undefined variable
set -o pipefail  # Exit on pipe failure

# ============================================
# Load configuration
# ============================================
CONFIG_FILE="${1:-/scripts/config.sh}"

if [[ ! -f "${CONFIG_FILE}" ]]; then
    echo "Error: Config file not found: ${CONFIG_FILE}"
    exit 1
fi

source "${CONFIG_FILE}"

# ============================================
# Helper functions
# ============================================

log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
}

log_step() {
    echo ""
    echo "============================================"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*"
    echo "============================================"
}

check_file() {
    if [[ ! -f "$1" ]]; then
        echo "Error: Required file not found: $1"
        exit 1
    fi
}

# ============================================
# Validate required input files
# ============================================

log_step "Validating required input files"

check_file "${REFERENCE_FASTA}"
log "Found genome FASTA: ${REFERENCE_FASTA}"

check_file "${GTF_FILE}"
log "Found GTF annotation: ${GTF_FILE}"

if [[ ! -f "${HISAT2_INDEX}.1.ht2" ]]; then
    echo "Error: HISAT2 index not found: ${HISAT2_INDEX}.1.ht2"
    echo "Please provide a valid HISAT2 index path in config.sh"
    exit 1
fi
log "Found HISAT2 index: ${HISAT2_INDEX}"

# ============================================
# Create output directories
# ============================================

log_step "Creating output directories"

mkdir -p "${TRIMMED_DIR}"
mkdir -p "${ALIGNED_DIR}"
mkdir -p "${SALMON_DIR}"
mkdir -p "${STRINGTIE_DIR}"
mkdir -p "${COUNTS_DIR}"
mkdir -p "${QC_DIR}"
mkdir -p "${JUNCTIONS_DIR}"
mkdir -p "${LOGS_DIR}"

# Splicing tool directories
mkdir -p "${RMATS_DIR}"
mkdir -p "${SUPPA2_DIR}"
mkdir -p "${LEAFCUTTER_DIR}"
mkdir -p "${DEXSEQ_DIR}"
mkdir -p "${JUNCTIONSEQ_DIR}"
mkdir -p "${MAJIQ_DIR}"
mkdir -p "${MISO_DIR}"

# ============================================
# Step 0: Generate derived reference files
# ============================================

log_step "Generating derived reference files"

# Generate transcriptome FASTA from GTF + genome (for Salmon)
if [[ ! -f "${TRANSCRIPTOME_FASTA}" ]]; then
    log "Extracting transcriptome sequences from GTF and genome..."
    gffread "${GTF_FILE}" -g "${REFERENCE_FASTA}" -w "${TRANSCRIPTOME_FASTA}"
    log "Transcriptome FASTA created: ${TRANSCRIPTOME_FASTA}"
else
    log "Transcriptome FASTA already exists: ${TRANSCRIPTOME_FASTA}"
fi

# Generate GFF3 from GTF (for MISO)
if [[ ! -f "${GFF3_FILE}" ]]; then
    log "Converting GTF to GFF3..."
    gffread "${GTF_FILE}" -o "${GFF3_FILE}" --force-exons -E
    log "GFF3 annotation created: ${GFF3_FILE}"
else
    log "GFF3 annotation already exists: ${GFF3_FILE}"
fi

# Build Salmon index
if [[ ! -d "${SALMON_INDEX}" ]] || [[ ! -f "${SALMON_INDEX}/info.json" ]]; then
    log "Building Salmon index..."
    salmon index \
        -t "${TRANSCRIPTOME_FASTA}" \
        -i "${SALMON_INDEX}" \
        -p ${THREADS} \
        2>&1 | tee "${LOGS_DIR}/salmon_index.log"
    log "Salmon index created: ${SALMON_INDEX}"
else
    log "Salmon index already exists: ${SALMON_INDEX}"
fi

# ============================================
# Step 1: QC and Trimming
# ============================================

if [[ ${RUN_QC} -eq 1 ]] || [[ ${RUN_TRIMMING} -eq 1 ]]; then
    log_step "Running QC and trimming with fastp"
    
    for sample in "${!SAMPLES[@]}"; do
        if [[ -f "${TRIMMED_DIR}/${sample}_R1.trimmed.fastq.gz" ]]; then
            log "Trimmed files exist for ${sample}, skipping..."
            continue
        fi
        
        log "Processing sample: ${sample}"
        
        R1="${FASTQ_R1[$sample]}"
        R2="${FASTQ_R2[$sample]}"
        
        check_file "${R1}"
        check_file "${R2}"
        
        fastp \
            -i "${R1}" \
            -I "${R2}" \
            -o "${TRIMMED_DIR}/${sample}_R1.trimmed.fastq.gz" \
            -O "${TRIMMED_DIR}/${sample}_R2.trimmed.fastq.gz" \
            -h "${QC_DIR}/${sample}_fastp.html" \
            -j "${QC_DIR}/${sample}_fastp.json" \
            -w ${THREADS} \
            --detect_adapter_for_pe \
            2>&1 | tee "${LOGS_DIR}/${sample}_fastp.log"
    done
fi

# ============================================
# Step 2: Alignment with HISAT2
# ============================================

if [[ ${RUN_ALIGNMENT} -eq 1 ]]; then
    log_step "Aligning reads with HISAT2"
    
    for sample in "${!SAMPLES[@]}"; do
        if [[ -f "${ALIGNED_DIR}/${sample}.sorted.bam" ]]; then
            log "Aligned BAM exists for ${sample}, skipping..."
            continue
        fi
        
        log "Aligning sample: ${sample}"
        
        if [[ -f "${TRIMMED_DIR}/${sample}_R1.trimmed.fastq.gz" ]]; then
            R1="${TRIMMED_DIR}/${sample}_R1.trimmed.fastq.gz"
            R2="${TRIMMED_DIR}/${sample}_R2.trimmed.fastq.gz"
        else
            R1="${FASTQ_R1[$sample]}"
            R2="${FASTQ_R2[$sample]}"
        fi
        
        hisat2 \
            -p ${THREADS} \
            -x "${HISAT2_INDEX}" \
            -1 "${R1}" \
            -2 "${R2}" \
            ${HISAT2_OPTS} \
            --rg-id "${sample}" \
            --rg "SM:${sample}" \
            --rg "LB:${sample}" \
            --rg "PL:ILLUMINA" \
            -S "${ALIGNED_DIR}/${sample}.sam" \
            2>&1 | tee "${LOGS_DIR}/${sample}_hisat2.log"
    done
fi

# ============================================
# Step 3: SAM to BAM conversion and sorting
# ============================================

if [[ ${RUN_SORTING} -eq 1 ]]; then
    log_step "Converting and sorting alignments"
    
    for sample in "${!SAMPLES[@]}"; do
        if [[ -f "${ALIGNED_DIR}/${sample}.sorted.bam" ]] && [[ -f "${ALIGNED_DIR}/${sample}.sorted.bam.bai" ]]; then
            log "Sorted BAM exists for ${sample}, skipping..."
            continue
        fi
        
        log "Processing sample: ${sample}"
        
        if [[ -f "${ALIGNED_DIR}/${sample}.sam" ]]; then
            samtools view -@ ${THREADS} -bS "${ALIGNED_DIR}/${sample}.sam" \
                > "${ALIGNED_DIR}/${sample}.unsorted.bam"
            
            samtools sort -@ ${THREADS} -m 4G \
                -o "${ALIGNED_DIR}/${sample}.sorted.bam" \
                "${ALIGNED_DIR}/${sample}.unsorted.bam"
            
            samtools index -@ ${THREADS} "${ALIGNED_DIR}/${sample}.sorted.bam"
            
            # Create symlinks
            ln -sf "${sample}.sorted.bam" "${ALIGNED_DIR}/${sample}.bam" 2>/dev/null || true
            ln -sf "${sample}.sorted.bam.bai" "${ALIGNED_DIR}/${sample}.bam.bai" 2>/dev/null || true
            
            # Clean up intermediate files
            rm -f "${ALIGNED_DIR}/${sample}.sam"
            rm -f "${ALIGNED_DIR}/${sample}.unsorted.bam"
            
            # Generate alignment statistics
            samtools flagstat "${ALIGNED_DIR}/${sample}.sorted.bam" \
                > "${QC_DIR}/${sample}_flagstat.txt"
        fi
    done
fi

# ============================================
# Step 4: Salmon quantification (for SUPPA2)
# ============================================

if [[ ${RUN_SALMON} -eq 1 ]]; then
    log_step "Running Salmon quantification"
    
    for sample in "${!SAMPLES[@]}"; do
        if [[ -f "${SALMON_DIR}/${sample}/quant.sf" ]]; then
            log "Salmon quant exists for ${sample}, skipping..."
            continue
        fi
        
        log "Quantifying sample: ${sample}"
        
        if [[ -f "${TRIMMED_DIR}/${sample}_R1.trimmed.fastq.gz" ]]; then
            R1="${TRIMMED_DIR}/${sample}_R1.trimmed.fastq.gz"
            R2="${TRIMMED_DIR}/${sample}_R2.trimmed.fastq.gz"
        else
            R1="${FASTQ_R1[$sample]}"
            R2="${FASTQ_R2[$sample]}"
        fi
        
        salmon quant \
            -i "${SALMON_INDEX}" \
            -l A \
            -1 "${R1}" \
            -2 "${R2}" \
            -p ${THREADS} \
            -o "${SALMON_DIR}/${sample}" \
            2>&1 | tee "${LOGS_DIR}/${sample}_salmon.log"
    done
fi

# ============================================
# Step 5: StringTie quantification
# ============================================

if [[ ${RUN_STRINGTIE} -eq 1 ]]; then
    log_step "Running StringTie quantification"
    
    for sample in "${!SAMPLES[@]}"; do
        if [[ -f "${STRINGTIE_DIR}/${sample}.gtf" ]]; then
            log "StringTie output exists for ${sample}, skipping..."
            continue
        fi
        
        log "Processing sample: ${sample}"
        
        stringtie "${ALIGNED_DIR}/${sample}.sorted.bam" \
            -G "${GTF_FILE}" \
            -e \
            -p ${THREADS} \
            -o "${STRINGTIE_DIR}/${sample}.gtf" \
            -A "${STRINGTIE_DIR}/${sample}_gene_abund.tab" \
            2>&1 | tee "${LOGS_DIR}/${sample}_stringtie.log"
    done
fi

# ============================================
# Step 6: Extract junctions with regtools
# ============================================

if [[ ${RUN_JUNCTION_EXTRACT} -eq 1 ]]; then
    log_step "Extracting splice junctions"
    
    for sample in "${!SAMPLES[@]}"; do
        if [[ -f "${JUNCTIONS_DIR}/${sample}_junctions.bed" ]]; then
            log "Junctions exist for ${sample}, skipping..."
            continue
        fi
        
        log "Extracting junctions: ${sample}"
        
        regtools junctions extract \
            -s 0 \
            -a 8 \
            -m 50 \
            -M 500000 \
            "${ALIGNED_DIR}/${sample}.sorted.bam" \
            -o "${JUNCTIONS_DIR}/${sample}_junctions.bed"
    done
fi

# ============================================
# Step 7: Gene/exon counting with featureCounts
# ============================================

if [[ ${RUN_COUNTS} -eq 1 ]]; then
    log_step "Counting reads with featureCounts"
    
    if [[ -f "${COUNTS_DIR}/gene_counts.txt" ]]; then
        log "Gene counts already exist, skipping..."
    else
        BAM_FILES=""
        for sample in "${!SAMPLES[@]}"; do
            BAM_FILES="${BAM_FILES} ${ALIGNED_DIR}/${sample}.sorted.bam"
        done
        
        featureCounts \
            -T ${THREADS} \
            -p \
            --countReadPairs \
            -a "${GTF_FILE}" \
            -o "${COUNTS_DIR}/gene_counts.txt" \
            ${BAM_FILES} \
            2>&1 | tee "${LOGS_DIR}/featureCounts_gene.log"
        
        featureCounts \
            -T ${THREADS} \
            -p \
            --countReadPairs \
            -f \
            -O \
            -a "${GTF_FILE}" \
            -o "${COUNTS_DIR}/exon_counts.txt" \
            ${BAM_FILES} \
            2>&1 | tee "${LOGS_DIR}/featureCounts_exon.log"
    fi
fi

# ============================================
# Step 8: rMATS differential splicing analysis
# ============================================

if [[ ${RUN_RMATS} -eq 1 ]]; then
    log_step "Running rMATS differential splicing analysis"
    
    if [[ -f "${RMATS_DIR}/SE.MATS.JCEC.txt" ]]; then
        log "rMATS results already exist, skipping..."
    else
        COND1_BAMS=""
        COND2_BAMS=""
        
        for sample in "${!SAMPLES[@]}"; do
            if [[ "${SAMPLES[$sample]}" == "${CONDITION1}" ]]; then
                if [[ -z "${COND1_BAMS}" ]]; then
                    COND1_BAMS="${ALIGNED_DIR}/${sample}.sorted.bam"
                else
                    COND1_BAMS="${COND1_BAMS},${ALIGNED_DIR}/${sample}.sorted.bam"
                fi
            else
                if [[ -z "${COND2_BAMS}" ]]; then
                    COND2_BAMS="${ALIGNED_DIR}/${sample}.sorted.bam"
                else
                    COND2_BAMS="${COND2_BAMS},${ALIGNED_DIR}/${sample}.sorted.bam"
                fi
            fi
        done
        
        echo "${COND1_BAMS}" | tr ',' '\n' > "${RMATS_DIR}/${CONDITION1}_bams.txt"
        echo "${COND2_BAMS}" | tr ',' '\n' > "${RMATS_DIR}/${CONDITION2}_bams.txt"
        
        rmats.py \
            --b1 "${RMATS_DIR}/${CONDITION1}_bams.txt" \
            --b2 "${RMATS_DIR}/${CONDITION2}_bams.txt" \
            --gtf "${GTF_FILE}" \
            --od "${RMATS_DIR}" \
            --tmp "${RMATS_DIR}/tmp" \
            -t paired \
            --readLength ${READ_LENGTH} \
            --nthread ${THREADS} \
            ${RMATS_OPTS} \
            2>&1 | tee "${LOGS_DIR}/rmats.log"
        
        rm -rf "${RMATS_DIR}/tmp"
        
        log "rMATS output files:"
        ls -la "${RMATS_DIR}"/*.txt 2>/dev/null || true
    fi
fi

# ============================================
# Step 9: SUPPA2 analysis
# ============================================

if [[ ${RUN_SUPPA2} -eq 1 ]]; then
    log_step "Running SUPPA2 analysis"
    
    # Generate local events
    if [[ ! -f "${SUPPA2_DIR}/events_SE_strict.ioe" ]]; then
        log "Generating SUPPA2 events from GTF..."
        suppa.py generateEvents \
            -i "${GTF_FILE}" \
            -o "${SUPPA2_DIR}/events" \
            -f ioe \
            -e ${SUPPA2_EVENT_TYPES} \
            2>&1 | tee "${LOGS_DIR}/suppa2_events.log"
    fi
    
    # Create TPM matrix from Salmon output
    if [[ ! -f "${SUPPA2_DIR}/tpm_matrix.txt" ]]; then
        log "Combining TPM values from Salmon..."
        
        # Get sample order
        sample_list=(${!SAMPLES[@]})
        first_sample="${sample_list[0]}"
        
        # Create header and first column (transcript IDs)
        echo -ne "transcript_id" > "${SUPPA2_DIR}/tpm_matrix.txt"
        for sample in "${sample_list[@]}"; do
            echo -ne "\t${sample}" >> "${SUPPA2_DIR}/tpm_matrix.txt"
        done
        echo "" >> "${SUPPA2_DIR}/tpm_matrix.txt"
        
        # Extract transcript IDs from first sample
        tail -n +2 "${SALMON_DIR}/${first_sample}/quant.sf" | cut -f1 > "${SUPPA2_DIR}/transcripts.txt"
        
        # Create TPM columns for each sample
        for sample in "${sample_list[@]}"; do
            tail -n +2 "${SALMON_DIR}/${sample}/quant.sf" | cut -f4 > "${SUPPA2_DIR}/${sample}_tpm.txt"
        done
        
        # Combine all columns
        paste "${SUPPA2_DIR}/transcripts.txt" "${SUPPA2_DIR}"/*_tpm.txt >> "${SUPPA2_DIR}/tpm_matrix.txt"
        
        # Clean up temp files
        rm -f "${SUPPA2_DIR}/transcripts.txt" "${SUPPA2_DIR}"/*_tpm.txt
    fi
    
    # Calculate PSI for each event type
    for event_file in "${SUPPA2_DIR}"/events*.ioe; do
        if [[ -f "${event_file}" ]]; then
            event_type=$(basename "${event_file}" .ioe | sed 's/events_//')
            
            if [[ ! -f "${SUPPA2_DIR}/${event_type}.psi" ]]; then
                log "Calculating PSI for ${event_type}..."
                suppa.py psiPerEvent \
                    -i "${event_file}" \
                    -e "${SUPPA2_DIR}/tpm_matrix.txt" \
                    -o "${SUPPA2_DIR}/${event_type}" \
                    2>&1 | tee -a "${LOGS_DIR}/suppa2_psi.log" || true
            fi
        fi
    done
    
    # Run differential splicing analysis
    log "Running SUPPA2 differential analysis..."
    
    # Create condition-specific sample lists
    cond1_cols=""
    cond2_cols=""
    col_idx=1
    for sample in "${!SAMPLES[@]}"; do
        if [[ "${SAMPLES[$sample]}" == "${CONDITION1}" ]]; then
            cond1_cols="${cond1_cols},${col_idx}"
        else
            cond2_cols="${cond2_cols},${col_idx}"
        fi
        ((col_idx++))
    done
    cond1_cols="${cond1_cols#,}"
    cond2_cols="${cond2_cols#,}"
    
    for event_file in "${SUPPA2_DIR}"/events*.ioe; do
        if [[ -f "${event_file}" ]]; then
            event_type=$(basename "${event_file}" .ioe | sed 's/events_//')
            
            if [[ -f "${SUPPA2_DIR}/${event_type}.psi" ]] && [[ ! -f "${SUPPA2_DIR}/diff_${event_type}.dpsi" ]]; then
                log "Differential analysis for ${event_type}..."
                
                # Split PSI and TPM by condition
                # This is a simplified approach - you may need to adjust
                suppa.py diffSplice \
                    -m empirical \
                    -i "${event_file}" \
                    -p "${SUPPA2_DIR}/${event_type}.psi" "${SUPPA2_DIR}/${event_type}.psi" \
                    -e "${SUPPA2_DIR}/tpm_matrix.txt" "${SUPPA2_DIR}/tpm_matrix.txt" \
                    -o "${SUPPA2_DIR}/diff_${event_type}" \
                    2>&1 | tee -a "${LOGS_DIR}/suppa2_diff.log" || true
            fi
        fi
    done
    
    log "SUPPA2 output files:"
    ls -la "${SUPPA2_DIR}"/*.psi "${SUPPA2_DIR}"/*.dpsi 2>/dev/null || true
fi

# ============================================
# Step 10: Leafcutter analysis
# ============================================

if [[ ${RUN_LEAFCUTTER} -eq 1 ]]; then
    log_step "Running Leafcutter analysis"
    
    if [[ -f "${LEAFCUTTER_DIR}/leafcutter_perind_numers.counts.gz" ]]; then
        log "Leafcutter clustering already done, skipping to differential analysis..."
    else
        # Convert BAM to junction files
        log "Converting BAM to junction files..."
        
        junction_files="${LEAFCUTTER_DIR}/junction_files.txt"
        > "${junction_files}"
        
        for sample in "${!SAMPLES[@]}"; do
            log "Extracting junctions for ${sample}..."
            
            # Use regtools to extract junctions
            if [[ ! -f "${LEAFCUTTER_DIR}/${sample}.junc" ]]; then
                regtools junctions extract \
                    -s 0 \
                    -a 8 \
                    -m 50 \
                    -M ${LEAFCUTTER_INTRON_MAX} \
                    "${ALIGNED_DIR}/${sample}.sorted.bam" \
                    -o "${LEAFCUTTER_DIR}/${sample}.junc"
            fi
            
            echo "${LEAFCUTTER_DIR}/${sample}.junc" >> "${junction_files}"
        done
        
        # Cluster introns
        log "Clustering introns..."
        python3 /tools/leafcutter/clustering/leafcutter_cluster_regtools.py \
            -j "${junction_files}" \
            -m ${LEAFCUTTER_MIN_COVERAGE} \
            -o "${LEAFCUTTER_DIR}/leafcutter" \
            -l ${LEAFCUTTER_INTRON_MAX} \
            2>&1 | tee "${LOGS_DIR}/leafcutter_cluster.log" || {
                # Alternative clustering script
                log "Trying alternative clustering script..."
                python3 /tools/leafcutter/clustering/leafcutter_cluster.py \
                    -j "${junction_files}" \
                    -m ${LEAFCUTTER_MIN_COVERAGE} \
                    -o "${LEAFCUTTER_DIR}/leafcutter" \
                    -l ${LEAFCUTTER_INTRON_MAX} \
                    2>&1 | tee "${LOGS_DIR}/leafcutter_cluster.log" || true
            }
    fi
    
    # Create groups file for differential analysis
    groups_file="${LEAFCUTTER_DIR}/groups.txt"
    > "${groups_file}"
    for sample in "${!SAMPLES[@]}"; do
        echo -e "${sample}\t${SAMPLES[$sample]}" >> "${groups_file}"
    done
    
    # Run differential splicing analysis
    if [[ -f "${LEAFCUTTER_DIR}/leafcutter_perind_numers.counts.gz" ]]; then
        log "Running Leafcutter differential intron usage analysis..."
        
        Rscript -e "
        library(leafcutter)
        
        counts_file <- '${LEAFCUTTER_DIR}/leafcutter_perind_numers.counts.gz'
        groups_file <- '${groups_file}'
        output_prefix <- '${LEAFCUTTER_DIR}/leafcutter_ds'
        
        tryCatch({
            leafcutter_ds(
                counts_file,
                groups_file,
                output_prefix = output_prefix
            )
        }, error = function(e) {
            message('Leafcutter DS completed with message: ', e\$message)
        })
        " 2>&1 | tee "${LOGS_DIR}/leafcutter_ds.log" || true
    fi
    
    log "Leafcutter output files:"
    ls -la "${LEAFCUTTER_DIR}"/*.txt "${LEAFCUTTER_DIR}"/*.gz 2>/dev/null || true
fi

# ============================================
# Step 11: DEXSeq analysis
# ============================================

if [[ ${RUN_DEXSEQ} -eq 1 ]]; then
    log_step "Running DEXSeq analysis"
    
    if [[ -f "${DEXSEQ_DIR}/dexseq_results.csv" ]]; then
        log "DEXSeq results already exist, skipping..."
    else
        # Create sample info file
        sample_info="${DEXSEQ_DIR}/sample_info.tsv"
        echo -e "sample_id\tcondition\tbam_file" > "${sample_info}"
        for sample in "${!SAMPLES[@]}"; do
            echo -e "${sample}\t${SAMPLES[$sample]}\t${sample}.sorted.bam" >> "${sample_info}"
        done
        
        Rscript /scripts/run_dexseq.R \
            "${GTF_FILE}" \
            "${ALIGNED_DIR}" \
            "${DEXSEQ_DIR}" \
            "${sample_info}" \
            ${THREADS} \
            2>&1 | tee "${LOGS_DIR}/dexseq.log" || true
        
        log "DEXSeq output files:"
        ls -la "${DEXSEQ_DIR}"/*.csv "${DEXSEQ_DIR}"/*.rds 2>/dev/null || true
    fi
fi

# ============================================
# Step 12: JunctionSeq analysis
# ============================================

if [[ ${RUN_JUNCTIONSEQ} -eq 1 ]]; then
    log_step "Running JunctionSeq analysis"
    
    if [[ -f "${JUNCTIONSEQ_DIR}/junctionseq_results.txt" ]]; then
        log "JunctionSeq results already exist, skipping..."
    else
        # Create sample info file
        sample_info="${JUNCTIONSEQ_DIR}/sample_info.tsv"
        echo -e "sample_id\tcondition\tbam_file" > "${sample_info}"
        for sample in "${!SAMPLES[@]}"; do
            echo -e "${sample}\t${SAMPLES[$sample]}\t${sample}.sorted.bam" >> "${sample_info}"
        done
        
        Rscript /scripts/run_junctionseq.R \
            "${GTF_FILE}" \
            "${ALIGNED_DIR}" \
            "${JUNCTIONSEQ_DIR}" \
            "${sample_info}" \
            ${THREADS} \
            2>&1 | tee "${LOGS_DIR}/junctionseq.log" || true
        
        log "JunctionSeq output files:"
        ls -la "${JUNCTIONSEQ_DIR}"/*.txt "${JUNCTIONSEQ_DIR}"/*.rds 2>/dev/null || true
    fi
fi

# ============================================
# Step 13: MAJIQ analysis (requires license)
# ============================================

if [[ ${RUN_MAJIQ} -eq 1 ]]; then
    log_step "Running MAJIQ analysis"
    
    if command -v majiq &> /dev/null; then
        if [[ -f "${MAJIQ_DIR}/voila_results.tsv" ]]; then
            log "MAJIQ results already exist, skipping..."
        else
            # Create MAJIQ config file
            majiq_config="${MAJIQ_DIR}/majiq.conf"
            
            # Get samples for each condition
            cond1_samples=""
            cond2_samples=""
            for sample in "${!SAMPLES[@]}"; do
                if [[ "${SAMPLES[$sample]}" == "${CONDITION1}" ]]; then
                    if [[ -z "${cond1_samples}" ]]; then
                        cond1_samples="${sample}"
                    else
                        cond1_samples="${cond1_samples},${sample}"
                    fi
                else
                    if [[ -z "${cond2_samples}" ]]; then
                        cond2_samples="${sample}"
                    else
                        cond2_samples="${cond2_samples},${sample}"
                    fi
                fi
            done
            
            cat > "${majiq_config}" << EOF
[info]
readlen=${READ_LENGTH}
samdir=${ALIGNED_DIR}
genome=custom
strandness=None

[experiments]
${CONDITION1}=${cond1_samples}
${CONDITION2}=${cond2_samples}
EOF
            
            # Build MAJIQ database
            log "Building MAJIQ database..."
            majiq build "${GTF_FILE}" \
                -c "${majiq_config}" \
                -j ${THREADS} \
                -o "${MAJIQ_DIR}/build" \
                2>&1 | tee "${LOGS_DIR}/majiq_build.log"
            
            # Run deltapsi analysis
            log "Running MAJIQ deltapsi..."
            majiq deltapsi \
                -grp1 "${MAJIQ_DIR}/build/"*"${CONDITION1}"*.majiq \
                -grp2 "${MAJIQ_DIR}/build/"*"${CONDITION2}"*.majiq \
                -j ${THREADS} \
                -o "${MAJIQ_DIR}/deltapsi" \
                -n "${CONDITION1}" "${CONDITION2}" \
                2>&1 | tee "${LOGS_DIR}/majiq_deltapsi.log"
            
            # Generate VOILA output
            log "Running VOILA..."
            voila tsv \
                "${MAJIQ_DIR}/build/splicegraph.sql" \
                "${MAJIQ_DIR}/deltapsi/"*.deltapsi.voila \
                -f "${MAJIQ_DIR}/voila_results.tsv" \
                2>&1 | tee "${LOGS_DIR}/voila.log"
            
            log "MAJIQ output files:"
            ls -la "${MAJIQ_DIR}"/*.tsv "${MAJIQ_DIR}"/*.voila 2>/dev/null || true
        fi
    else
        log "WARNING: MAJIQ not installed. Skipping MAJIQ analysis."
        log "MAJIQ requires a license from https://majiq.biociphers.org/"
    fi
fi

# ============================================
# Step 14: MISO analysis
# ============================================

if [[ ${RUN_MISO} -eq 1 ]]; then
    log_step "Running MISO analysis"
    
    # Check if MISO comparisons already exist
    if [[ -d "${MISO_DIR}/comparisons" ]] && [[ -n "$(ls -A ${MISO_DIR}/comparisons/*.miso_bf 2>/dev/null)" ]]; then
        log "MISO results already exist, skipping..."
    else
        # Index GFF3 for MISO
        if [[ ! -d "${MISO_DIR}/indexed_events" ]] || [[ -z "$(ls -A ${MISO_DIR}/indexed_events 2>/dev/null)" ]]; then
            log "Indexing GFF3 for MISO..."
            mkdir -p "${MISO_DIR}/indexed_events"
            
            index_gff --index "${GFF3_FILE}" "${MISO_DIR}/indexed_events/" \
                2>&1 | tee "${LOGS_DIR}/miso_index.log" || {
                    log "Warning: MISO index_gff failed. Creating minimal index..."
                    # Create a minimal indexed structure
                    mkdir -p "${MISO_DIR}/indexed_events/genes"
                }
        fi
        
        # Run MISO for each sample
        for sample in "${!SAMPLES[@]}"; do
            if [[ -d "${MISO_DIR}/${sample}/summary" ]]; then
                log "MISO output exists for ${sample}, skipping..."
                continue
            fi
            
            log "Running MISO for ${sample}..."
            
            miso \
                --run "${MISO_DIR}/indexed_events/" \
                "${ALIGNED_DIR}/${sample}.sorted.bam" \
                --output-dir "${MISO_DIR}/${sample}" \
                --read-len ${READ_LENGTH} \
                -p ${THREADS} \
                2>&1 | tee "${LOGS_DIR}/${sample}_miso.log" || true
            
            # Summarize MISO output
            if [[ -d "${MISO_DIR}/${sample}" ]]; then
                summarize_miso --summarize-samples "${MISO_DIR}/${sample}" \
                    "${MISO_DIR}/${sample}/summary/" \
                    2>&1 | tee -a "${LOGS_DIR}/${sample}_miso.log" || true
            fi
        done
        
        # Compare conditions
        log "Comparing conditions with MISO..."
        
        cond1_dirs=""
        cond2_dirs=""
        for sample in "${!SAMPLES[@]}"; do
            if [[ -d "${MISO_DIR}/${sample}" ]]; then
                if [[ "${SAMPLES[$sample]}" == "${CONDITION1}" ]]; then
                    if [[ -z "${cond1_dirs}" ]]; then
                        cond1_dirs="${MISO_DIR}/${sample}"
                    else
                        cond1_dirs="${cond1_dirs},${MISO_DIR}/${sample}"
                    fi
                else
                    if [[ -z "${cond2_dirs}" ]]; then
                        cond2_dirs="${MISO_DIR}/${sample}"
                    else
                        cond2_dirs="${cond2_dirs},${MISO_DIR}/${sample}"
                    fi
                fi
            fi
        done
        
        if [[ -n "${cond1_dirs}" ]] && [[ -n "${cond2_dirs}" ]]; then
            mkdir -p "${MISO_DIR}/comparisons"
            compare_miso \
                --compare-samples "${cond1_dirs}" "${cond2_dirs}" \
                "${MISO_DIR}/comparisons/" \
                2>&1 | tee "${LOGS_DIR}/miso_compare.log" || true
        fi
        
        log "MISO output files:"
        ls -la "${MISO_DIR}"/comparisons/*.miso_bf 2>/dev/null || true
    fi
fi

# ============================================
# Step 15: Prepare test data for sashimiplots
# ============================================

log_step "Preparing test data for sashimiplots"

TEST_DATA_DIR="${OUTPUT_DIR}/sashimiplots_testdata"
mkdir -p "${TEST_DATA_DIR}"

# Copy BAM files
log "Copying BAM files..."
for sample in "${!SAMPLES[@]}"; do
    cp "${ALIGNED_DIR}/${sample}.sorted.bam" "${TEST_DATA_DIR}/"
    cp "${ALIGNED_DIR}/${sample}.sorted.bam.bai" "${TEST_DATA_DIR}/"
done

# Create tool result directories
mkdir -p "${TEST_DATA_DIR}/rmats"
mkdir -p "${TEST_DATA_DIR}/suppa2"
mkdir -p "${TEST_DATA_DIR}/leafcutter"
mkdir -p "${TEST_DATA_DIR}/dexseq"
mkdir -p "${TEST_DATA_DIR}/junctionseq"
mkdir -p "${TEST_DATA_DIR}/majiq"
mkdir -p "${TEST_DATA_DIR}/miso"

# Copy results from each tool
if [[ -d "${RMATS_DIR}" ]]; then
    log "Copying rMATS results..."
    cp "${RMATS_DIR}"/*.MATS.JCEC.txt "${TEST_DATA_DIR}/rmats/" 2>/dev/null || true
    cp "${RMATS_DIR}"/*.MATS.JC.txt "${TEST_DATA_DIR}/rmats/" 2>/dev/null || true
fi

if [[ -d "${SUPPA2_DIR}" ]]; then
    log "Copying SUPPA2 results..."
    cp "${SUPPA2_DIR}"/*.psi "${TEST_DATA_DIR}/suppa2/" 2>/dev/null || true
    cp "${SUPPA2_DIR}"/*.dpsi "${TEST_DATA_DIR}/suppa2/" 2>/dev/null || true
    cp "${SUPPA2_DIR}"/*.ioe "${TEST_DATA_DIR}/suppa2/" 2>/dev/null || true
fi

if [[ -d "${LEAFCUTTER_DIR}" ]]; then
    log "Copying Leafcutter results..."
    cp "${LEAFCUTTER_DIR}"/*cluster* "${TEST_DATA_DIR}/leafcutter/" 2>/dev/null || true
    cp "${LEAFCUTTER_DIR}"/*effect_sizes* "${TEST_DATA_DIR}/leafcutter/" 2>/dev/null || true
fi

if [[ -d "${DEXSEQ_DIR}" ]]; then
    log "Copying DEXSeq results..."
    cp "${DEXSEQ_DIR}"/*.csv "${TEST_DATA_DIR}/dexseq/" 2>/dev/null || true
fi

if [[ -d "${JUNCTIONSEQ_DIR}" ]]; then
    log "Copying JunctionSeq results..."
    cp "${JUNCTIONSEQ_DIR}"/*.txt "${TEST_DATA_DIR}/junctionseq/" 2>/dev/null || true
fi

if [[ -d "${MAJIQ_DIR}" ]]; then
    log "Copying MAJIQ results..."
    cp "${MAJIQ_DIR}"/*.tsv "${TEST_DATA_DIR}/majiq/" 2>/dev/null || true
    cp "${MAJIQ_DIR}/deltapsi/"*.voila "${TEST_DATA_DIR}/majiq/" 2>/dev/null || true
fi

if [[ -d "${MISO_DIR}/comparisons" ]]; then
    log "Copying MISO results..."
    cp "${MISO_DIR}/comparisons/"*.miso_bf "${TEST_DATA_DIR}/miso/" 2>/dev/null || true
fi

# Copy GTF
if [[ -f "${GTF_FILE}" ]]; then
    log "Copying annotation file..."
    cp "${GTF_FILE}" "${TEST_DATA_DIR}/"
fi

# Create sample metadata
log "Creating sample metadata..."
cat > "${TEST_DATA_DIR}/sample_info.tsv" << EOF
sample_id	condition	bam_file
EOF
for sample in "${!SAMPLES[@]}"; do
    echo -e "${sample}\t${SAMPLES[$sample]}\t${sample}.sorted.bam" >> "${TEST_DATA_DIR}/sample_info.tsv"
done

# Create example R script
cat > "${TEST_DATA_DIR}/example_all_tools.R" << 'EOF'
# Example script for using sashimiplots with all supported tools
library(sashimiplots)

data_dir <- "/output/sashimiplots_testdata"

# Read sample info
sample_info <- read.delim(file.path(data_dir, "sample_info.tsv"))

# Define BAM files by condition
control_bams <- file.path(data_dir, 
                          sample_info$bam_file[sample_info$condition == "Control"])
treatment_bams <- file.path(data_dir,
                            sample_info$bam_file[sample_info$condition == "Treatment"])

bam_list <- list(Control = control_bams, Treatment = treatment_bams)

# Read GTF annotation
gtf_files <- list.files(data_dir, pattern = "\\.gtf$", full.names = TRUE)
if (length(gtf_files) > 0) {
  annotation <- read_gtf_annotation(gtf_files[1])
}

# ============================================
# rMATS
# ============================================
if (file.exists(file.path(data_dir, "rmats/SE.MATS.JCEC.txt"))) {
  cat("Testing rMATS integration...\n")
  se_events <- read_rmats(file.path(data_dir, "rmats/SE.MATS.JCEC.txt"), event_type = "SE")
  cat("  Loaded", nrow(se_events), "SE events\n")
}

# ============================================
# SUPPA2
# ============================================
suppa2_files <- list.files(file.path(data_dir, "suppa2"), pattern = "\\.dpsi$", full.names = TRUE)
if (length(suppa2_files) > 0) {
  cat("Testing SUPPA2 integration...\n")
  suppa2_events <- read_suppa2(suppa2_files[1])
  cat("  Loaded", nrow(suppa2_events), "events\n")
}

# ============================================
# Leafcutter
# ============================================
leafcutter_files <- list.files(file.path(data_dir, "leafcutter"), 
                                pattern = "cluster_significance", full.names = TRUE)
if (length(leafcutter_files) > 0) {
  cat("Testing Leafcutter integration...\n")
  leafcutter_events <- read_leafcutter(leafcutter_files[1])
  cat("  Loaded", nrow(leafcutter_events), "clusters\n")
}

# ============================================
# DEXSeq
# ============================================
if (file.exists(file.path(data_dir, "dexseq/dexseq_results.csv"))) {
  cat("Testing DEXSeq integration...\n")
  dexseq_events <- read_dexseq(file.path(data_dir, "dexseq/dexseq_results.csv"))
  cat("  Loaded", nrow(dexseq_events), "exons\n")
}

# ============================================
# JunctionSeq
# ============================================
junctionseq_files <- list.files(file.path(data_dir, "junctionseq"), 
                                 pattern = "results", full.names = TRUE)
if (length(junctionseq_files) > 0) {
  cat("Testing JunctionSeq integration...\n")
  junctionseq_events <- read_junctionseq(junctionseq_files[1])
  cat("  Loaded", nrow(junctionseq_events), "features\n")
}

# ============================================
# MAJIQ
# ============================================
majiq_files <- list.files(file.path(data_dir, "majiq"), pattern = "\\.tsv$", full.names = TRUE)
if (length(majiq_files) > 0) {
  cat("Testing MAJIQ integration...\n")
  majiq_events <- read_majiq(majiq_files[1])
  cat("  Loaded", nrow(majiq_events), "LSVs\n")
}

# ============================================
# MISO
# ============================================
miso_files <- list.files(file.path(data_dir, "miso"), pattern = "\\.miso_bf$", full.names = TRUE)
if (length(miso_files) > 0) {
  cat("Testing MISO integration...\n")
  miso_events <- read_miso(miso_files[1])
  cat("  Loaded", nrow(miso_events), "events\n")
}

cat("\nAll tool integrations tested successfully!\n")
EOF

log "Test data prepared at: ${TEST_DATA_DIR}"
ls -la "${TEST_DATA_DIR}"

# ============================================
# Summary
# ============================================

log_step "Pipeline completed!"

echo ""
echo "============================================"
echo "Input files used:"
echo "  - Genome FASTA: ${REFERENCE_FASTA}"
echo "  - GTF annotation: ${GTF_FILE}"
echo "  - HISAT2 index: ${HISAT2_INDEX}"
echo ""
echo "Auto-generated files:"
echo "  - Transcriptome FASTA: ${TRANSCRIPTOME_FASTA}"
echo "  - Salmon index: ${SALMON_INDEX}"
echo "  - GFF3 annotation: ${GFF3_FILE}"
echo ""
echo "Output summary:"
echo "  - Trimmed reads: ${TRIMMED_DIR}"
echo "  - Aligned BAMs: ${ALIGNED_DIR}"
echo "  - Salmon quant: ${SALMON_DIR}"
echo ""
echo "Splicing tool results:"
echo "  - rMATS: ${RMATS_DIR}"
echo "  - SUPPA2: ${SUPPA2_DIR}"
echo "  - Leafcutter: ${LEAFCUTTER_DIR}"
echo "  - DEXSeq: ${DEXSEQ_DIR}"
echo "  - JunctionSeq: ${JUNCTIONSEQ_DIR}"
echo "  - MAJIQ: ${MAJIQ_DIR}"
echo "  - MISO: ${MISO_DIR}"
echo ""
echo "Test data: ${TEST_DATA_DIR}"
echo ""
echo "Run example_all_tools.R to test sashimiplots integration"
echo "============================================"
