FROM mambaorg/micromamba:0.15.3

# apt installs
USER root
# nuke cache dirs before installing pkgs; tip from Dirk E fixes broken img
RUN rm -f /var/lib/dpkg/available && rm -rf  /var/cache/apt/*
RUN apt-get update && \
  apt-get -y upgrade && \
  apt-get install -y --no-install-recommends \
  procps \
  python2.7 \
  build-essential curl git python3-pip procps \ 
  manpages-dev g++ gcc gfortran make autoconf automake libtool \
  zlib1g-dev liblzma-dev libbz2-dev lbzip2 libgsl-dev \
  libblas-dev libx11-dev \
  libreadline-dev libxt-dev libpcre2-dev libcurl4-openssl-dev \
  && rm -rf /var/lib/apt/lists/*

# install mamba env
USER micromamba
COPY --chown=micromamba:micromamba environment.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

ENV MAMBA_ROOT_PREFIX="/opt/conda"
WORKDIR "$MAMBA_ROOT_PREFIX"

RUN which pip
# Add additional software using Conda env:
RUN pip install cellSNP && pip install vireoSNP
    
# test main python libraries can be loaded:
RUN python -c 'import sys;print(sys.version_info);import celltypist; import cellSNP; import vireoSNP; import scanpy; import click; import pandas; import plotnine; import matplotlib'

## check software versions:
RUN celltypist --help  >> /opt/conda/conda_software_versions.txt 2>&1 || true
RUN cellSNP  >> /opt/conda/conda_software_versions.txt 2>&1 || true
RUN cellsnp-lite -V >> /opt/conda/conda_software_versions.txt 2>&1 || true
RUN vireo  >> /opt/conda/conda_software_versions.txt 2>&1 || true
RUN which python >> /opt/conda/conda_software_versions.txt
RUN python --version >> /opt/conda/conda_software_versions.txt
RUN tabix --version  >> /opt/conda/conda_software_versions.txt
RUN bcftools --version  >> /opt/conda/conda_software_versions.txt
RUN samtools --version  >> /opt/conda/conda_software_versions.txt
RUN cat /opt/conda/conda_software_versions.txt

USER micromamba
WORKDIR /tmp
ENTRYPOINT ["/bin/entrypoint.sh"]
CMD ["/bin/bash"]