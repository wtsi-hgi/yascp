FROM python:3.10.2-buster
# run critical system updates ... one could also use the artful-curl pack

RUN apt-get update \
    && apt-get install -y \
    apt-utils \
    build-essential \
    git \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-gnutls-dev \
    libssl-dev \
    libperl-dev \
    libgsl0-dev \
    && apt-get purge \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}

WORKDIR /opt
RUN wget https://github.com/samtools/htslib/releases/download/1.15.1/htslib-1.15.1.tar.bz2 \
  && tar jxvf htslib-1.15.1.tar.bz2 \
  && cd htslib-1.15.1 \
  && ./configure \
  && make \
  && make install

WORKDIR /opt
RUN wget https://github.com/samtools/bcftools/releases/download/1.15.1/bcftools-1.15.1.tar.bz2 \
  && tar jxvf bcftools-1.15.1.tar.bz2 \
  && cd bcftools-1.15.1 \
  && ./configure \
  && make \
  && make install

WORKDIR /opt
FROM r-base:4.2.0

WORKDIR /
