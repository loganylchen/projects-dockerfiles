FROM mambaorg/micromamba:latest
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV PATH=/opt/conda/bin:$PATH
ARG VERSION=0.69.9

RUN micromamba install -y -c conda-forge -c bioconda circos=${VERSION} git\
    && micromamba clean --all --yes