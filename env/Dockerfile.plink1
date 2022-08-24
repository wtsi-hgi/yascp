#FROM ubuntu:18.04
FROM python:3.10.2-buster

# Installation instructions: https://www.cog-genomics.org/plink
# Environment variables, plink 1.9
ENV PLINK_VERSION       20220402
ENV PLINK_HOME          /usr/local/plink
ENV PATH                $PLINK_HOME:$PATH

RUN apt-get update \
    && apt-get install -y \
    apt-utils \
    unzip \
    wget \
    && apt-get purge \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* && \
    apt-get autoremove -y && \
    rm -rf /var/lib/{apt,dpkg,cache,log}

RUN wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_${PLINK_VERSION}.zip && \
  unzip plink_linux_x86_64_${PLINK_VERSION}.zip -d $PLINK_HOME && \
  rm plink_linux_x86_64_${PLINK_VERSION}.zip 

# Set the default action to print plink's options 
CMD ["plink"]
