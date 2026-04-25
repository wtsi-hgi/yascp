# Use the latest available Miniforge image
FROM condaforge/miniforge3:latest

# Set environment variables
ENV PATH="/opt/conda/bin:${PATH}" \
    DEBIAN_FRONTEND=noninteractive \
    TZ=Etc/UTC \
    NB_USER=container_user

# Update and install base system dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    wget \
    unzip \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libgsl-dev \
    libboost-all-dev \
    dirmngr \
    apt-transport-https \
    ca-certificates \
    software-properties-common \
    gnupg2 \
    curl \
    git \
    g++ \
    gcc \
    gfortran \
    make \
    autoconf \
    automake \
    libtool \
    zlib1g-dev \
    liblzma-dev \
    libbz2-dev \
    lbzip2 \
    libblas-dev \
    libx11-dev \
    libreadline-dev \
    libxt-dev \
    libpcre2-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*
RUN apt-get update && \
    apt-get install --yes -y --no-install-recommends dirmngr apt-transport-https ca-certificates software-properties-common gnupg2 \
    build-essential curl git \ 
    g++ gcc gfortran make autoconf automake libtool \
    zlib1g-dev liblzma-dev libbz2-dev lbzip2 libgsl-dev \
    libblas-dev libx11-dev libncurses-dev wget \
    libreadline-dev libxt-dev libpcre2-dev libcurl4-openssl-dev && \
    rm -rf /var/lib/apt/lists/*

# # Install BCFtools
RUN wget https://github.com/samtools/bcftools/releases/download/1.21/bcftools-1.21.tar.bz2 \
  && tar jxvf bcftools-1.21.tar.bz2 \
  && cd bcftools-1.21 \
  && ./configure \
  && make \
  && make install

# Install HTSlib
RUN mkdir -p /opt/htslib && \
    cd /opt/htslib && \
    wget -q https://github.com/samtools/htslib/releases/download/1.3.1/htslib-1.3.1.tar.bz2 && \
    tar -xjf htslib-1.3.1.tar.bz2 && \
    cd htslib-1.3.1 && make && \
    rm -f /opt/htslib/htslib-1.3.1.tar.bz2

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    libcurl4-openssl-dev \
    libssl-dev \
    libxml2-dev \
    libgit2-dev \
    libgsl-dev \
    libboost-all-dev \
    libx11-dev \
    libxt-dev \
    libreadline-dev \
    libpcre2-dev && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2 \
  && tar jxvf samtools-1.18.tar.bz2 \
  && cd samtools-1.18 \
  && ./configure \
  && make \
  && make install

# Copy the environment.yml file into the container
COPY environment.yml /tmp/environment.yml

# Create the Conda environment from the YAML file
RUN conda config --set always_yes yes --set changeps1 no && conda env create -f /tmp/environment.yml && \
    conda clean -afy

RUN apt-get install -y libxml2-dev libcurl4-openssl-dev libssl-dev make openssl libxml2
RUN apt install dirmngr gnupg apt-transport-https ca-certificates software-properties-common
RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9
RUN wget -qO- https://cloud.r-project.org/bin/linux/ubuntu/marutter_pubkey.asc | gpg --dearmor -o /usr/share/keyrings/cran-archive-keyring.gpg
RUN echo "deb [signed-by=/usr/share/keyrings/cran-archive-keyring.gpg] https://cloud.r-project.org/bin/linux/ubuntu $(lsb_release -cs)-cran40/" > /etc/apt/sources.list.d/cran-r.list

RUN apt update --yes; apt upgrade --yes
RUN apt-get update && apt-get install -y \
    libpng-dev libtiff5-dev libjpeg-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfreetype6-dev \
    pkg-config \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Activate the environment in every shell session
ENV PATH="/opt/conda/envs/yascp-qc-env/bin:$PATH"

# Set default Conda environment for Python
RUN echo "conda activate yascp-qc-env" >> ~/.bashrc

SHELL ["conda", "run", "-n", "yascp-qc-env", "/bin/bash", "-c"]

RUN R --version
# InstallR and necessary R packages
RUN which R
RUN apt-get update && apt-get install -y libblas-dev gfortran && apt-get install -y liblapack-dev libblas-dev gfortran libv8-dev wget libhdf5-dev libcairo2-dev libz-dev
RUN apt install -y libglpk-dev libgmp-dev libopenblas-dev liblapack-dev

# Step 1: Install core scientific packages
RUN conda install numpy pytables hdf5 python-igraph leidenalg numba llvmlite
RUN conda install -c conda-forge tensorflow
RUN conda update plotnine
# Step 2: Install PyTorch and related libraries
RUN conda install -c pytorch -c nvidia pytorch torchvision torchaudio pytorch-cuda=11.7
RUN conda clean --all --yes && pip cache purge
# Step 3: Install bioinformatics tools
RUN conda install -c conda-forge hdf5
RUN conda install anndata scanpy matplotlib
# Step 3: Use Pip for remaining Python packages
RUN pip install doubletdetection cellex scvi-tools cellSNP vireoSNP xtarfile
# Step 4: Clean up Conda and Pip cache
RUN mkdir -p /tmp/numba_cache && chmod -R 777 /tmp/numba_cache

RUN python -c "import numpy; print(numpy.__version__)"
RUN python -c "import tables; print(tables.__version__)"
RUN conda install -c conda-forge zlib libcurl libxml2 libpng
RUN conda install -c conda-forge r-curl r-httpuv r-png r-data.table r-httr r-shiny r-igraph r-reticulate r-plotly r-miniui r-leiden r-fastdummies
RUN conda install -c conda-forge r-ragg
RUN conda install bioconductor-fgsea
RUN conda install -c conda-forge bioconductor-dose
RUN conda install bioconductor-TFBSTools
RUN conda install r-cairo
RUN conda install -c conda-forge libcurl 
RUN conda install -c pytorch -c nvidia pytorch torchvision torchaudio pytorch-cuda=11.7
RUN conda install -c conda-forge jax[cuda12] jaxlib
RUN pip install --upgrade jax jaxlib==0.4.29+cuda12.cudnn91 -f https://storage.googleapis.com/jax-releases/jax_cuda_releases.html
RUN conda install -c conda-forge scvi-tools

RUN conda install  \
    r-ggpubr r-nloptr \
    r-units \
    r-odbc \
    r-git2r \
    r-ragg \
    r-rmariadb \
    r-lme4 \
    r-sf \
    r-epir \
    r-mediation \
    r-spdep \
    r-spatialreg \
    r-car \
    r-aer \
    r-drc \
    r-systemfit
RUN conda install -c bioconda -c conda-forge celltypist


#  Works till here ok. 
RUN R -e "install.packages('tidyverse', dependencies = TRUE, repos = 'https://cloud.r-project.org')"
RUN R --no-echo -e 'install.packages(c(  "readr", "cpp11", "broom", "dbplyr", "ggplot2", "haven", "lubridate", "modelr", "ragg", "readr", "readxl", "tidyr","igraph","remotes", "optparse","data.table", "ggplot2", "devtools", "BiocManager"), repos=c("http://cran.us.r-project.org","https://cran.ma.imperial.ac.uk/", "https://www.stats.bris.ac.uk/R/"), dependencies = TRUE)'
RUN R --no-echo -e 'install.packages(c( "joineRML"), dependencies = TRUE,repos=c("https://cloud.r-project.org","https://cran.ma.imperial.ac.uk/", "https://www.stats.bris.ac.uk/R/"))'
RUN R --no-echo -e 'install.packages(c( "gtable", "ggplot2", "cowplot", "RcppArmadillo", "RcppProgress","harmony"), dependencies = TRUE,repos=c("https://cloud.r-project.org","https://cran.ma.imperial.ac.uk/", "https://www.stats.bris.ac.uk/R/"))'
RUN R --no-echo -e "remotes::install_github(c('YuLab-SMU/enrichplot'), dependencies = TRUE)"
RUN R --no-echo -e 'install.packages(c( "argparse"), repos = "https://cloud.r-project.org")'
RUN R --no-echo -e 'BiocManager::install(c("clusterProfiler", "org.Hs.eg.db" ,"DeconRNASeq","limma","scDblFinder", "clusterProfiler", "hopach", "mygene", "tidyr", "R.utils", "foreach", "doParallel", "stringr"), dependencies = TRUE)'
RUN R --no-echo -e "remotes::install_github(c('EDePasquale/DoubletDecon'), dependencies = TRUE)"
RUN R --no-echo -e 'install.packages("pkgdown", repos = "https://cloud.r-project.org", dependencies = TRUE)'
RUN R --no-echo -e 'install.packages("dsb", repos = "https://cloud.r-project.org", dependencies = TRUE)'
RUN R -e "install.packages(c('ggrastr'), repos = 'https://cloud.r-project.org', dependencies = TRUE)"
RUN R --quiet -e 'install.packages(c("Matrix","igraph","leiden"), repos = "https://cloud.r-project.org", dependencies = TRUE)'
RUN R --quiet -e "remotes::install_github('satijalab/seurat', ref = 'seurat5', quiet = TRUE)"
RUN R --no-echo -e "remotes::install_github(c('mojaveazure/seurat-disk'), dependencies = TRUE)"
RUN R --no-echo -e "remotes::install_github('satijalab/azimuth', ref = 'master')"
RUN R --no-echo -e "remotes::install_github('EDePasquale/DoubletDecon', ref = 'master')"
RUN R --no-echo -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder', ref = 'master')"
RUN R --no-echo -e "remotes::install_github(c('mojaveazure/seurat-disk'), dependencies = TRUE)"

# RUN which conda
RUN mkdir /opt/PBMC_reference
ADD https://zenodo.org/record/4546839/files/idx.annoy https://zenodo.org/record/4546839/files/ref.Rds /opt/PBMC_reference/
RUN chmod 755 -R /opt/PBMC_reference
# Verify that harfbuzz and fribidi `.pc` files are present
RUN find /usr -name "harfbuzz.pc" || echo "harfbuzz.pc not found"
RUN find /usr -name "fribidi.pc" || echo "fribidi.pc not found"
RUN find /usr -name "freetype2.pc" || echo "freetype2.pc not found"
RUN find /usr -name "libpng.pc" || echo "libpng.pc not found"
# Install tidyverse
RUN R --no-echo -e 'BiocManager::install(c("scDblFinder"), dependencies = TRUE)'
RUN R --no-echo -e 'BiocManager::install(c("scds"), dependencies = TRUE)'
RUN R --no-echo -e "remotes::install_github(c('joseah/HierscPred'), dependencies = TRUE)"
RUN R --version
RUN conda update scanpy


# Verify installation
RUN pip install scikeras
RUN python -c "from scikeras.wrappers import KerasClassifier, KerasRegressor"

