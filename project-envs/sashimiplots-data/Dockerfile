# RNA-seq Analysis Environment for sashimiplots Test Data Generation
# Includes all splicing analysis tools supported by sashimiplots:
# - rMATS, SUPPA2, Leafcutter, DEXSeq, JunctionSeq, MAJIQ, MISO

FROM ubuntu:22.04

LABEL maintainer="Logan Chen <loganylchen@gmail.com>"
LABEL description="RNA-seq analysis pipeline for generating sashimiplots test data"

# Prevent interactive prompts
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ="Asia/Hong_Kong"

# ARG
ARG RMATS_VERSION=4.3.0

# PATH
ENV PATH=/opt/bin:/opt/rmats/:${PATH}

# Install system dependencies
RUN apt-get update && apt-get install -y \
    wget \
    curl \
    git \
    build-essential \
    cmake \
    autoconf \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libncurses5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    python3 \
    python3-pip \
    python3-dev \
    python3-venv \
    default-jre \
    unzip \
    gzip \
    pigz \
    parallel \
    bedtools \
    gcc build-essential   python3-pip \
    ca-certificates \
    cmake \
    curl \
    cython3 \
    g++ \
    gfortran \
    git \
    libblas-dev \
    libgsl-dev \
    liblapack-dev \
    make \
    python-is-python3 \
    python3-dev \
    r-base \
    r-cran-doparallel \
    r-cran-dosnow \
    r-cran-foreach \
    r-cran-getopt \
    r-cran-ggplot2 \
    r-cran-iterators \
    r-cran-mixtools \
    r-cran-nloptr \
    r-cran-rcpp \
    zlib1g-dev \
    && update-alternatives --install /usr/bin/python python /usr/bin/python3.10 1 \
    && update-alternatives --install /usr/bin/pip pip /usr/bin/pip3 1 \
    && rm -rf /var/lib/apt/lists/*

# Create working directories
RUN mkdir -p /tools /data /output /reference /scripts

# ============================================
# Install gffread (for extracting transcriptome from GTF + genome)
# ============================================
WORKDIR /tools
RUN git clone https://github.com/gpertea/gffread.git && \
    cd gffread && \
    make release && \
    cp gffread /usr/local/bin/



# Install R packages
RUN R -e "install.packages(c('BiocManager', 'devtools', 'data.table', 'dplyr', 'ggplot2', 'optparse'), repos='https://cloud.r-project.org/')"

# Install Bioconductor packages for DEXSeq and JunctionSeq
RUN R -e "BiocManager::install(c('DEXSeq', 'JunctionSeq', 'GenomicRanges', 'GenomicFeatures', 'Rsamtools', 'rtracklayer', 'GenomicAlignments'), ask=FALSE, update=FALSE)"

# ============================================
# Install HISAT2
# ============================================
WORKDIR /tools
RUN wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download -O hisat2-2.2.1-Linux_x86_64.zip && \
    unzip hisat2-2.2.1-Linux_x86_64.zip && \
    rm hisat2-2.2.1-Linux_x86_64.zip
ENV PATH="/tools/hisat2-2.2.1:${PATH}"

# ============================================
# Install samtools
# ============================================
RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 && \
    tar -xjf samtools-1.18.tar.bz2 && \
    cd samtools-1.18 && \
    ./configure --prefix=/tools/samtools && \
    make -j4 && \
    make install && \
    cd .. && rm -rf samtools-1.18 samtools-1.18.tar.bz2
ENV PATH="/tools/samtools/bin:${PATH}"

# ============================================
# Install HTSlib (for bgzip, tabix)
# ============================================
RUN wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2 && \
    tar -xjf htslib-1.18.tar.bz2 && \
    cd htslib-1.18 && \
    ./configure --prefix=/tools/htslib && \
    make -j4 && \
    make install && \
    cd .. && rm -rf htslib-1.18 htslib-1.18.tar.bz2
ENV PATH="/tools/htslib/bin:${PATH}"

# ============================================
# Install StringTie (for transcript assembly and TPM)
# ============================================
RUN wget http://ccb.jhu.edu/software/stringtie/dl/stringtie-2.2.1.Linux_x86_64.tar.gz && \
    tar -xzf stringtie-2.2.1.Linux_x86_64.tar.gz && \
    rm stringtie-2.2.1.Linux_x86_64.tar.gz
ENV PATH="/tools/stringtie-2.2.1.Linux_x86_64:${PATH}"

# ============================================
# Install Python packages
# ============================================
RUN pip3 install --no-cache-dir \
    numpy \
    scipy \
    pysam \
    Cython \
    pandas \
    matplotlib \
    seaborn \
    scikit-learn \
    statsmodels

# ============================================
# Install rMATS (differential splicing analysis)
# ============================================
RUN cd /tmp/ && \
    mkdir rmats_build \
    && cd rmats_build \
    && git clone https://github.com/Xinglab/rmats-turbo.git \
    && cd rmats-turbo \
    && git checkout v${RMATS_VERSION} \
    && echo '' > setup_environment.sh \
    && ./build_rmats \
    && mkdir -p rmats \
    && cd rmats \
    && cp /tmp/rmats_build/rmats-turbo/rmats.py ./ \
    && cp /tmp/rmats_build/rmats-turbo/cp_with_prefix.py ./ \
    && cp /tmp/rmats_build/rmats-turbo/*.so ./ \
    && mkdir rMATS_C \
    && cp /tmp/rmats_build/rmats-turbo/rMATS_C/rMATSexe ./rMATS_C \
    && mkdir rMATS_P \
    && cp /tmp/rmats_build/rmats-turbo/rMATS_P/*.py ./rMATS_P \
    && mkdir rMATS_R \
    && cp /tmp/rmats_build/rmats-turbo/rMATS_R/*.R ./rMATS_R \
    && cp -r ../rmats /opt \
    && chmod -R +x /opt

# ============================================
# Install SUPPA2 (alternative splicing analysis)
# ============================================
RUN pip3 install --no-cache-dir SUPPA==2.3

# ============================================
# Install Leafcutter
# ============================================
WORKDIR /tools
RUN git clone https://github.com/davidaknowles/leafcutter.git && \
    cd leafcutter && \
    pip3 install -r requirements.txt || true
ENV PATH="/tools/leafcutter/scripts:/tools/leafcutter/clustering:${PATH}"

# Install Leafcutter R package
RUN R -e "devtools::install_github('davidaknowles/leafcutter/leafcutter')"

# ============================================
# Install MISO (Mixture of Isoforms)
# ============================================
RUN pip3 install --no-cache-dir misopy

# ============================================
# Install MAJIQ (Manual installation required - license needed)
# Note: MAJIQ requires a license from https://majiq.biociphers.org/
# This section creates placeholder for manual installation
# ============================================
RUN mkdir -p /tools/majiq && \
    echo "MAJIQ requires a license. Download from https://majiq.biociphers.org/" > /tools/majiq/README.txt
# If you have MAJIQ wheel file, uncomment and modify:
# COPY majiq-*.whl /tools/majiq/
# RUN pip3 install /tools/majiq/majiq-*.whl

# ============================================
# Install regtools (junction extraction)
# ============================================
WORKDIR /tools
RUN git clone https://github.com/griffithlab/regtools.git && \
    cd regtools && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make && \
    cp regtools /usr/local/bin/
ENV PATH="/tools/regtools/build:${PATH}"

# ============================================
# Install fastp (read QC and trimming)
# ============================================
RUN wget http://opengene.org/fastp/fastp && \
    chmod a+x fastp && \
    mv fastp /usr/local/bin/

# ============================================
# Install featureCounts (read counting)
# ============================================
RUN cd /tools/ && \
    wget https://github.com/ShiLab-Bioinformatics/subread/releases/download/2.0.2/subread-2.0.2-Linux-x86_64.tar.gz && \
    tar -xzf subread-2.0.2-Linux-x86_64.tar.gz && \
    rm subread-2.0.2-Linux-x86_64.tar.gz
ENV PATH="/tools/subread-2.0.6-Linux-x86_64/bin:${PATH}"

# ============================================
# Install Salmon (for TPM quantification - needed by SUPPA2)
# ============================================
RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz && \
    tar -xzf salmon-1.10.0_linux_x86_64.tar.gz && \
    rm salmon-1.10.0_linux_x86_64.tar.gz
ENV PATH="/tools/salmon-latest_linux_x86_64/bin:${PATH}"

# ============================================
# Create helper scripts directory
# ============================================
WORKDIR /scripts

# Copy analysis scripts
COPY run_analysis.sh /scripts/
COPY config.sh /scripts/
COPY run_dexseq.R /scripts/
COPY run_junctionseq.R /scripts/
RUN chmod +x /scripts/*.sh

# Set working directory
WORKDIR /data

# Default command
CMD ["/bin/bash"]
