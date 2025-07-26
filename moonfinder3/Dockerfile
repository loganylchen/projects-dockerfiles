FROM btrspg/vscode-base:0.0.6

ADD install_packages.R /tmp/
ADD bioconductor_packages.txt /tmp/
ADD devtools_packages.txt /tmp/
ADD CRAN_packages.txt /tmp/
ADD requirements.txt /tmp/


RUN Rscript /tmp/install_packages.R  && \
    pip install -r /tmp/requirements.txt && \
    useradd -m -u 501 logan && \
    apt autoremove && \
    apt clean && \
    apt autoclean && \
    rm -rf /tmp/*