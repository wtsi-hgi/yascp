FROM ubuntu:18.04
RUN cat /etc/os-release
ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && \
    apt-get install --yes -y --no-install-recommends dirmngr apt-transport-https ca-certificates software-properties-common gnupg2 \
    build-essential curl git \ 
    g++ gcc gfortran make autoconf automake libtool \
    zlib1g-dev liblzma-dev libbz2-dev lbzip2 libgsl-dev \
    libblas-dev libx11-dev \
    libreadline-dev libxt-dev libpcre2-dev libcurl4-openssl-dev && \
    rm -rf /var/lib/apt/lists/*
RUN apt install dirmngr gnupg apt-transport-https ca-certificates software-properties-common
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran40/'
RUN apt update --yes; apt upgrade --yes
RUN apt-get install

RUN apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev make openssl libxml2
RUN apt install r-base --yes
RUN apt install build-essential
RUN apt install apt-transport-https software-properties-common
RUN R --version
RUN which R

RUN R -e ".libPaths(); "
RUN R -e "install.packages(c('ggrastr'), dependencies = TRUE)"

# # RUN R --no-echo -e ".libPaths('/usr/local/lib/R/site-library'); install.packages(c('ragg'), dependencies = TRUE)"
RUN apt install libcurl4-gnutls-dev --yes
# install.packages(c("devtools", "roxygen2", "pkgdown", "covr")); devtools::install_dev_deps()
# RUN apt install r-cran-ragg r-cran-xml2
RUN apt install libtiff-dev  --yes
RUN apt-get install libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev build-essential libcurl4-openssl-dev libxml2-dev libssl-dev libfontconfig1-dev --yes
RUN R --no-echo -e 'install.packages("pkgdown")'
# RUN R --no-echo -e 'install.packages(c("devtools", "roxygen2","covr")); devtools::install_dev_deps()'
RUN R --no-echo -e "install.packages(c('tidyverse'))"
RUN R --no-echo -e "install.packages(c('remotes'))"
RUN R --no-echo -e 'library(remotes); remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)'
# RUN Rscript --no-echo -e ".libPaths('/usr/local/lib/R/site-library'); install.packages(c('lifecycle'), repos = 'https://cran.r-project.org', dependencies = TRUE)"
# RUN R --no-echo -e ".libPaths('/usr/local/lib/R/site-library'); install.packages(c('tidyverse'), dependencies = TRUE)"
RUN apt-get install -y libv8-dev
RUN apt-get install -y wget

RUN apt-get install libhdf5-dev --yes
RUN R --no-echo -e "install.packages(c('ggplot2', 'hdf5r'), repo='https://cloud.r-project.org')"
RUN R --no-echo -e 'install.packages("hdf5r")'
RUN R --no-echo -e "remotes::install_github('satijalab/azimuth', ref = 'master')"
RUN R --no-echo -e "remotes::install_github(c('mojaveazure/seurat-disk'), dependencies = TRUE)"
RUN ls -lh /usr/local/lib/R/site-library
# RUN which conda
RUN mkdir /opt/PBMC_reference
ADD https://zenodo.org/record/4546839/files/idx.annoy https://zenodo.org/record/4546839/files/ref.Rds /opt/PBMC_reference/
RUN chmod 755 -R /opt/PBMC_reference
RUN R --no-echo -e 'install.packages(c("devtools"))'
RUN apt-get install libcairo2-dev  --yes
# RUN R --no-echo -e "devtools::install_github('VPetukhov/ggrastr', build_vignettes = TRUE)"
RUN R -e "install.packages(c('ggrastr'), dependencies = TRUE)"
RUN R --no-echo -e "install.packages('BiocManager'); BiocManager::install('limma')"
RUN R --no-echo -e "install.packages('dsb')"