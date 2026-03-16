# RNA-seq Analysis Pipeline for sashimiplots Test Data

This directory contains a Docker-based RNA-seq analysis pipeline for generating test data compatible with **all** splicing tools supported by sashimiplots.

## Supported Tools

| Tool | Description | Output Format |
|------|-------------|---------------|
| **rMATS** | Replicate Multivariate Analysis of Transcript Splicing | `*.MATS.JCEC.txt` |
| **SUPPA2** | Super-fast PSI calculation and differential analysis | `*.dpsi`, `*.psi` |
| **Leafcutter** | Intron cluster differential splicing | `*cluster_significance.txt` |
| **DEXSeq** | Differential exon usage | `dexseq_results.csv` |
| **JunctionSeq** | Differential junction usage | `junctionseq_results.txt` |
| **MAJIQ** | Local Splicing Variation detection | `*.voila.tsv` |
| **MISO** | Mixture of Isoforms | `*.miso_bf` |

## Quick Start

### 1. Build Docker Image

```bash
cd scripts/analysis
docker build -t rnaseq-sashimi-analysis .
```

### 2. Prepare Your Data

```
/path/to/your/data/
├── fastq/
│   ├── Sample1_R1.fastq.gz    # Condition 1, Replicate 1
│   ├── Sample1_R2.fastq.gz
│   ├── Sample2_R1.fastq.gz    # Condition 1, Replicate 2
│   ├── Sample2_R2.fastq.gz
│   ├── Sample3_R1.fastq.gz    # Condition 1, Replicate 3
│   ├── Sample3_R2.fastq.gz
│   ├── Sample4_R1.fastq.gz    # Condition 2, Replicate 1
│   ├── Sample4_R2.fastq.gz
│   ├── Sample5_R1.fastq.gz    # Condition 2, Replicate 2
│   ├── Sample5_R2.fastq.gz
│   ├── Sample6_R1.fastq.gz    # Condition 2, Replicate 3
│   └── Sample6_R2.fastq.gz
├── reference/
│   ├── genome.fa              # Reference genome
│   ├── transcriptome.fa       # Transcriptome (for Salmon)
│   └── annotation.gtf         # GTF annotation
└── output/                    # Created by pipeline
```

### 3. Modify Configuration

Edit `config.sh`:

```bash
# Sample names and conditions
SAMPLES["Sample1"]="Control"
SAMPLES["Sample2"]="Control"
SAMPLES["Sample3"]="Control"
SAMPLES["Sample4"]="Treatment"
SAMPLES["Sample5"]="Treatment"
SAMPLES["Sample6"]="Treatment"

# FASTQ paths
FASTQ_R1["Sample1"]="/data/fastq/Sample1_R1.fastq.gz"
FASTQ_R2["Sample1"]="/data/fastq/Sample1_R2.fastq.gz"
# ... etc

# Enable/disable tools
RUN_RMATS=1
RUN_SUPPA2=1
RUN_LEAFCUTTER=1
RUN_DEXSEQ=1
RUN_JUNCTIONSEQ=1
RUN_MAJIQ=0      # Requires license
RUN_MISO=1
```

### 4. Run Analysis

```bash
docker run -it --rm \
    -v /path/to/fastq:/data/fastq \
    -v /path/to/reference:/reference \
    -v /path/to/output:/output \
    -v $(pwd)/config.sh:/scripts/config.sh \
    rnaseq-sashimi-analysis \
    /scripts/run_analysis.sh
```

## Output Structure

```
/output/
├── aligned/                    # Sorted BAM files
├── salmon/                     # Salmon quantification
├── stringtie/                  # StringTie quantification
├── counts/                     # featureCounts gene/exon counts
├── rmats/                      # rMATS results
│   ├── SE.MATS.JCEC.txt       # Skipped exon events
│   ├── A3SS.MATS.JCEC.txt     # Alt 3' splice site
│   ├── A5SS.MATS.JCEC.txt     # Alt 5' splice site
│   ├── MXE.MATS.JCEC.txt      # Mutually exclusive exons
│   └── RI.MATS.JCEC.txt       # Retained introns
├── suppa2/                     # SUPPA2 results
│   ├── events_SE.ioe          # Event definitions
│   ├── SE.psi                 # PSI values
│   └── diff_SE.dpsi           # Differential PSI
├── leafcutter/                 # Leafcutter results
│   └── leafcutter_cluster_significance.txt
├── dexseq/                     # DEXSeq results
│   ├── dexseq_results.csv
│   └── dexseq_significant.csv
├── junctionseq/                # JunctionSeq results
│   └── junctionseq_results.txt
├── majiq/                      # MAJIQ results (if licensed)
│   └── voila_results.tsv
├── miso/                       # MISO results
│   └── comparisons/*.miso_bf
└── sashimiplots_testdata/      # Ready-to-use test data
    ├── *.sorted.bam
    ├── rmats/
    ├── suppa2/
    ├── leafcutter/
    ├── dexseq/
    ├── junctionseq/
    ├── majiq/
    ├── miso/
    ├── sample_info.tsv
    └── example_all_tools.R
```

## Using with sashimiplots

```r
library(sashimiplots)

# Read different tool outputs
se_events <- read_rmats("rmats/SE.MATS.JCEC.txt", event_type = "SE")
suppa_events <- read_suppa2("suppa2/diff_SE.dpsi")
leafcutter_clusters <- read_leafcutter("leafcutter/cluster_significance.txt")
dexseq_results <- read_dexseq("dexseq/dexseq_results.csv")
junctionseq_results <- read_junctionseq("junctionseq/junctionseq_results.txt")
majiq_lsvs <- read_majiq("majiq/voila_results.tsv")
miso_events <- read_miso("miso/comparisons/comparison.miso_bf")

# Plot sashimi for any significant event
region <- get_rmats_region(se_events[1, ], padding = 500)
plot_sashimi(
  bam_files = list(Control = control_bams, Treatment = treatment_bams),
  region = region,
  annotation = annotation
)
```

## Tool-Specific Notes

### MAJIQ
MAJIQ requires a license from https://majiq.biociphers.org/. Set `RUN_MAJIQ=0` in config.sh or install the wheel file manually in the Docker image.

### MISO
MISO requires properly formatted GFF3 annotations. The pipeline attempts automatic conversion from GTF, but manual adjustment may be needed for complex annotations.

### Leafcutter
Leafcutter uses regtools for junction extraction. Ensure BAM files are properly sorted and indexed.

## Configuration Options

| Parameter | Description | Default |
|-----------|-------------|---------|
| `THREADS` | CPU threads | 8 |
| `READ_LENGTH` | Read length | 150 |
| `CONDITION1` | First condition name | "Control" |
| `CONDITION2` | Second condition name | "Treatment" |
| `LEAFCUTTER_MIN_COVERAGE` | Min junction reads for Leafcutter | 30 |
| `SUPPA2_EVENT_TYPES` | SUPPA2 event types | "SE SS MX RI FL" |
| `FDR_CUTOFF` | FDR threshold for DEXSeq/JunctionSeq | 0.05 |

## Pipeline Control Flags

```bash
# Core pipeline
RUN_QC=1
RUN_TRIMMING=1
RUN_ALIGNMENT=1
RUN_SORTING=1
RUN_SALMON=1
RUN_STRINGTIE=1
RUN_COUNTS=1
RUN_JUNCTION_EXTRACT=1

# Splicing tools
RUN_RMATS=1
RUN_SUPPA2=1
RUN_LEAFCUTTER=1
RUN_DEXSEQ=1
RUN_JUNCTIONSEQ=1
RUN_MAJIQ=0          # Disabled by default (requires license)
RUN_MISO=1
```

## Troubleshooting

### Memory Issues
- Reduce `THREADS` in config.sh
- Increase Docker memory limit: `docker run --memory=32g ...`

### Tool-specific Errors
Check individual log files in `/output/logs/`:
- `rmats.log`
- `suppa2_*.log`
- `leafcutter_*.log`
- `dexseq.log`
- `junctionseq.log`
- `miso_*.log`

### Missing Dependencies
If a tool fails, you can rebuild the Docker image or run tools individually inside the container.

## License

MIT License - see main sashimiplots package.
# sashimiplots-data-docker
