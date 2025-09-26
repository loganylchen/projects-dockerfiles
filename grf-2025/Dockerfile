FROM mambaorg/micromamba:latest
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:$PATH
ARG VERSION=0.69.9

ARG USER_UID=1777
ARG USER_GID=1777
ARG USER_NAME=loganvscode


RUN micromamba install -y -c conda-forge -c bioconda circos=${VERSION} git stringtie r-irkernel nanocount r-languageserver \
    && git config --global user.name "loganylchen" \
    && git config --global user.email "yuelong.chen.btr@gmail.com"  \
    && micromamba clean --all --yes



ADD install_packages.R /tmp/
ADD bioconductor_packages.txt /tmp/
ADD devtools_packages.txt /tmp/
ADD CRAN_packages.txt /tmp/
ADD requirements.txt /tmp/


RUN Rscript /tmp/install_packages.R  && \
    pip install -r /tmp/requirements.txt 

ADD this.R /tmp/
RUN Rscript /tmp/this.R
RUN groupadd --gid $USER_GID loganvscode && \
    useradd -u $USER_UID -g $USER_NAME -s /bin/bash -m $USER_NAME  && \
    apt autoremove && \
    apt clean && \
    apt autoclean && \
    rm -rf /tmp/*

USER $USER_NAME 
