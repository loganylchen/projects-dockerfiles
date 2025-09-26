FROM mambaorg/micromamba:1.3.1
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:$PATH
ENV MAMBA_ROOT_PREFIX=/opt/conda
ENV ENV_NAME=GRF
ARG VERSION=0.69.9

ARG USER_UID=1777
ARG USER_GID=1777
ARG USER_NAME=loganvscode

USER root

RUN micromamba create -y -n ${ENV_NAME} git r-base=4.4.3 python=3.10  &&\
    micromamba install -n ${ENV_NAME} -y -c conda-forge -c bioconda circos=${VERSION} \
    && git config --global user.name "loganylchen" \
    && git config --global user.email "yuelong.chen.btr@gmail.com"  \
    && micromamba clean --all --yes



ADD install_tools.bash /tmp/
ADD CRAN_packages.bash /tmp/
ADD devtools_packages.R /tmp/
ADD bioconductor_packages.bash /tmp/
ADD requirements.txt /tmp/



RUN bash /tmp/install_tools.bash && \
    bash /tmp/CRAN_packages.bash && \
    bash /tmp/bioconductor_packages.bash && \
    micromamba run -n ${ENV_NAME} Rscript /tmp/devtools_packages.R && \
    micromamba run -n ${ENV_NAME} pip install -r /tmp/requirements.txt 




ADD this.R /tmp/
RUN Rscript /tmp/this.R



RUN groupadd --gid $USER_GID loganvscode && \
    useradd -u $USER_UID -g $USER_NAME -s /bin/bash -m $USER_NAME  && \
    apt autoremove && \
    apt clean && \
    apt autoclean && \
    rm -rf /tmp/*

USER $USER_NAME 

CMD ["micromamba", "run", "-n", "${ENV_NAME}", "bash"]

