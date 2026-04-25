FROM us.gcr.io/broad-dsde-methods/cellbender:latest
# ARG conda_env variable must match conda env name in environment.yml:
ARG conda_env=conda

LABEL authors="Matiss Ozols" \
  maintainer="Matiss Ozols <mo11@sanger.ac.uk>" \
  description="Docker image for WSTI-HGI scRNA cellbender pipeline"

# nuke cache dirs before installing pkgs; tip from Dirk E fixes broken img
RUN rm -f /var/lib/dpkg/available && rm -rf  /var/cache/apt/*
RUN apt-get update && \
    apt-get install --yes --no-install-recommends \
        wget \
        bzip2 \
        ca-certificates \
        curl \
        git \
        zip \
        unzip \
        procps && \
    apt-get install --yes \
        libpng-dev \
        libcurl4-gnutls-dev \
        libssl-dev \
        libxml2-dev \
        libgit2-dev \
        zlib1g-dev \
        build-essential \
  	procps && \ 
     apt-get purge && \
     apt-get clean && \
     rm -rf /var/lib/apt/lists/*

ENV CONDA_DEFAULT_ENV $conda_env
RUN which cellbender
RUN git clone https://github.com/broadinstitute/CellBender.git
RUN ls
RUN pip install -e CellBender
RUN pip install --upgrade cellbender==0.3.2
RUN which cellbender

RUN /opt/conda/bin/cellbender -h
RUN pip install plotnine
RUN python -c 'import sys;print(sys.version_info);import click; import pandas; import plotnine; import matplotlib'
RUN pip show cellbender

CMD /bin/sh