FROM continuumio/miniconda3
# the folllowing ARG conda_env variable must match the conda env name defined in environment.yml:
ARG conda_env=bff

LABEL authors="M.Ozols" \
  maintainer="M.Ozols <mo11@sanger.ak>" \
  description="Docker image for samtools, freemuxlet, bedtools, freebayes"

# nuke cache dirs before installing pkgs; tip from Dirk E fixes broken img
RUN rm -f /var/lib/dpkg/available && rm -rf  /var/cache/apt/*
RUN apt-get update && \
  apt-get -y upgrade && \
  apt-get install -y --no-install-recommends \
  r-base \
  python2.7 \
  build-essential curl git procps \ 
  g++ gcc gfortran make autoconf automake libtool \
  zlib1g-dev liblzma-dev libbz2-dev lbzip2 libgsl-dev \
  libblas-dev libx11-dev \
  libreadline-dev libxt-dev libpcre2-dev libcurl4-openssl-dev \
  && rm -rf /var/lib/apt/lists/*
RUN  apt update

# update conda && install Conda env:
RUN conda update -n base -c defaults conda
ADD environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml
ENV CONDA_DEFAULT_ENV $conda_env
ENV PATH /opt/conda/envs/$conda_env/bin:$PATH

CMD /bin/sh
