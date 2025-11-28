FROM mambaorg/micromamba:latest

ARG USER_UID=1777
ARG USER_GID=1777
ARG USER_NAME=loganvscode
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:$PATH
ENV ANTHROPIC_BASE_URL=https://api.moonshot.cn/anthropic
ENV ANTHROPIC_MODEL=kimi-k2-thinking-turbo
ENV ANTHROPIC_DEFAULT_OPUS_MODEL=kimi-k2-thinking-turbo
ENV ANTHROPIC_DEFAULT_SONNET_MODEL=kimi-k2-thinking-turbo
ENV ANTHROPIC_DEFAULT_HAIKU_MODEL=kimi-k2-thinking-turbo
ENV CLAUDE_CODE_SUBAGENT_MODEL=kimi-k2-thinking-turbo 

RUN micromamba install -y \
    conda-forge::r-ggalt \
    conda-forge::r-forestplot \
    conda-forge::git \
    anaconda::jupyter \
    conda-forge::r-vip \
    bioconda::bioconductor-annotationdbi \
    bioconda::bioconductor-biobase \
    bioconda::bioconductor-biocstyle \
    bioconda::bioconductor-biomart \
    bioconda::bioconductor-bsgenome.hsapiens.ucsc.hg38 \
    bioconda::bioconductor-clusterprofiler \
    bioconda::bioconductor-complexheatmap \
    bioconda::bioconductor-deseq2 \
    bioconda::bioconductor-dittoseq \
    bioconda::bioconductor-edger \
    bioconda::bioconductor-enhancedvolcano \
    bioconda::bioconductor-genefilter \
    bioconda::bioconductor-geoquery \
    bioconda::bioconductor-pd.hg.u133.plus.2 \
    bioconda::bioconductor-hgu133plus2.db \
    bioconda::bioconductor-affy \
    conda-forge::r-ggextra \
    conda-forge::libopenblas=0.3.6 \
    bioconda::bioconductor-glmgampoi \
    bioconda::bioconductor-gosemsim \
    bioconda::bioconductor-oligo \
    bioconda::bioconductor-gseabase \
    bioconda::bioconductor-gsva \
    bioconda::bioconductor-impute \
    bioconda::bioconductor-limma \
    bioconda::bioconductor-mast \
    bioconda::bioconductor-org.hs.eg.db \
    bioconda::bioconductor-pcatools \
    bioconda::bioconductor-rcy3 \
    bioconda::bioconductor-scater \ 
    conda-forge::r-biocmanager \
    conda-forge::r-circlize \
    conda-forge::r-cluster \
    conda-forge::r-corrplot\
    conda-forge::r-devtools \
    conda-forge::r-dplyr \
    conda-forge::r-extrafont \
    conda-forge::r-ggextra \
    conda-forge::r-ggplot2 \
    conda-forge::r-ggplotify \
    conda-forge::r-ggpubr \
    conda-forge::r-ggrepel \
    conda-forge::r-ggsci \
    conda-forge::r-ggvenndiagram \
    conda-forge::r-glmnet \
    conda-forge::r-vegan \
    conda-forge::r-igraph \
    conda-forge::r-irkernel \
    conda-forge::r-knitr \
    conda-forge::r-languageserver \
    conda-forge::r-pheatmap \
    conda-forge::r-proc \
    conda-forge::r-progress \
    conda-forge::r-prroc \
    conda-forge::r-purrr \
    conda-forge::r-rcolorbrewer \
    bioconda::bioconductor-rcy3 \
    conda-forge::r-readr \
    conda-forge::r-rlang \
    conda-forge::r-rmarkdown \
    conda-forge::r-stringr \
    conda-forge::r-testthat \
    conda-forge::r-tibble \
    conda-forge::r-tidyr \
    conda-forge::r-tidyverse \
    conda-forge::r-usethis \
    && micromamba clean --all --yes

ADD ./install_packages.R /tmp/
ADD ./devtools_packages.txt /tmp/


USER root
RUN micromamba run -n base Rscript /tmp/install_packages.R &&\
    git config --global user.name "loganylchen" && \
    git config --global user.email "yuelong.chen.btr@gmail.com" && \
    groupadd --gid $USER_GID loganvscode && \
    useradd -u $USER_UID -g $USER_NAME -s /bin/bash -m $USER_NAME  && \
    chmod 777 -R /opt/conda/lib/R/library && \
    apt autoremove && \
    apt clean && \
    apt autoclean && \
    rm -rf /tmp/*
USER $USER_NAME 

ADD install_claude_code.bash /tmp/ 
RUN bash /tmp/install_claude_code.bash && rm -rf /tmp/*
