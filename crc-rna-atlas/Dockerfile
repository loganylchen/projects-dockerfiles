FROM mambaorg/micromamba:latest

ARG USER_UID=1777
ARG USER_GID=1777
ARG USER_NAME=loganvscode
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:$PATH


RUN micromamba install -y \
    conda-forge::r-ggplot2 \
    conda-forge::r-biocmanager \
    conda-forge::r-dplyr \
    conda-forge::r-ggpubr \
    conda-forge::r-ggsci \
    conda-forge::r-tidyr \
    conda-forge::r-tibble \
    conda-forge::r-devtools \
    conda-forge::r-readr \
    conda-forge::r-tidyverse \
    conda-forge::r-ggrepel \
    conda-forge::r-ggplotify \
    conda-forge::r-igraph \
    conda-forge::r-cluster \
    conda-forge::r-ggplot2 \
    conda-forge::r-dplyr \
    conda-forge::r-tidyr \
    conda-forge::r-tibble \
    conda-forge::r-readr \
    conda-forge::r-stringr \
    conda-forge::r-purrr \
    conda-forge::r-rlang \
    bioconda::bioconductor-biobase \
    bioconda::bioconductor-org.hs.eg.db \
    bioconda::bioconductor-annotationdbi \
    conda-forge::r-testthat \
    conda-forge::r-knitr \
    conda-forge::r-rmarkdown \
    bioconda::bioconductor-biocstyle \
    conda-forge::r-pheatmap \
    conda-forge::r-rcolorbrewer \
    conda-forge::r-circlize \
    conda-forge::r-languageserver \
    conda-forge::r-irkernel \
    conda-forge::r-usethis \
    conda-forge::r-progress \
    conda-forge::r-glmnet \
    conda-forge::r-proc \
    conda-forge::r-ggpubr \
    conda-forge::r-ggsci \
    conda-forge::r-extrafont \
    conda-forge::r-ggvenndiagram \
    conda-forge::r-prroc \
    conda-forge::r-corrplot\
    conda-forge::r-ggextra \
    bioconda::bioconductor-gseabase \
    bioconda::bioconductor-gsva \
    bioconda::bioconductor-enhancedvolcano \
    bioconda::bioconductor-pcatools \
    bioconda::bioconductor-deseq2 \
    bioconda::bioconductor-edger \
    bioconda::bioconductor-limma \
    bioconda::bioconductor-pcatools \
    bioconda::bioconductor-mast \
    bioconda::bioconductor-org.hs.eg.db \
    bioconda::bioconductor-dittoseq \
    bioconda::bioconductor-scater \
    bioconda::bioconductor-glmgampoi \
    bioconda::bioconductor-clusterprofiler \
    bioconda::bioconductor-glmgampoi \
    bioconda::bioconductor-complexheatmap \
    bioconda::bioconductor-annotationdbi \
    bioconda::bioconductor-biomart \
    bioconda::bioconductor-geoquery \
    bioconda::bioconductor-gosemsim \
    bioconda::bioconductor-gsva \
    bioconda::bioconductor-genefilter \
    bioconda::bioconductor-rcy3 \
    bioconda::bioconductor-impute \
    && micromamba clean --all --yes

ADD ./install_packages.R /tmp/
ADD ./devtools_packages.txt /tmp/


USER root
RUN micromamba run -n base Rscript /tmp/install_packages.R &&\
    git config --global user.name "loganylchen" && \
    git config --global user.email "yuelong.chen.btr@gmail.com" && \
    groupadd --gid $USER_GID loganvscode && \
    useradd -u $USER_UID -g $USER_NAME -s /bin/bash -m $USER_NAME  && \
    apt autoremove && \
    apt clean && \
    apt autoclean && \
    rm -rf /tmp/*
USER $USER_NAME 
