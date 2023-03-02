FROM wtsihgi/nf_scrna_qc:6bb6af5
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
RUN chown -R ${NB_USER}:${NB_USER} ${HOME}
USER ${NB_USER}
SHELL ["conda", "run", "-n", "cenv", "/bin/bash", "-c"]

# RUN which pip
RUN pip install scanpy==1.9.2
############################ Setup: R packages #################################
ENV CONDA_DEFAULT_ENV cenv
ENV PATH ${HOME}/conda/envs/cenv/bin:$PATH
USER root

USER ${NB_USER}
################################################################################
