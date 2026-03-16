FROM btrspg/vscode-base:0.0.6

ARG USER_UID=1777
ARG USER_GID=1777
ARG USER_NAME=loganvscode

ADD install_packages.R /tmp/
ADD bioconductor_packages.txt /tmp/
ADD devtools_packages.txt /tmp/
ADD CRAN_packages.txt /tmp/
ADD requirements.txt /tmp/


RUN Rscript /tmp/install_packages.R  && \
    pip install -r /tmp/requirements.txt 
RUN git config --global user.name "loganylchen" && \
    git config --global user.email "yuelong.chen.btr@gmail.com" && \
    groupadd --gid $USER_GID loganvscode && \
    useradd -u $USER_UID -g $USER_NAME -s /bin/bash -m $USER_NAME  && \
    apt autoremove && \
    apt clean && \
    apt autoclean && \
    rm -rf /tmp/*
USER $USER_NAME 
