FROM mercury/nf_qc_scrna:v2
# For this version we have altered the scanpy version to 1.9.2 as there is an issue that Onek1k dataset expereineces in pca step. 
USER root
############################ Setup: OS #########################################
# Set up locales properly
# RUN apt-get update && \
#     apt-get install --yes --no-install-recommends locales && \
#     apt-get purge && \
#     apt-get clean && \
#     rm -rf /var/lib/apt/lists/*
ENV NB_USER container_user
ENV HOME /home/${NB_USER}

# USER container_user
# RUN chown -R ${NB_USER}:${NB_USER} ${HOME}
USER ${NB_USER}

RUN conda config --add channels defaults & conda config --add channels bioconda & conda config --add channels conda-forge & conda config --set channel_priority strict
RUN ls
# RUN conda install openssl
# RUN conda install  -c bioconda samtools bcftools plink 
RUN pip install doubletdetection
RUN pip install  matplotlib xtarfile

USER ${NB_USER}
################################################################################
# import doubletdetection
# import tarfile
# import matplotlib
# matplotlib.use('PDF')
# import matplotlib.pyplot as plt
# import os
# import argparse
# import sys
# import pandas as pd