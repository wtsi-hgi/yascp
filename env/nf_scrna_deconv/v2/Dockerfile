FROM mercury/scrna_deconvolution:62bd56a
# For this version we have altered the scanpy version to 1.9.2 as there is an issue that Onek1k dataset expereineces in pca step. 
USER micromamba
############################ Setup: OS #########################################
# Set up locales properly
# RUN apt-get update && \
#     apt-get install --yes --no-install-recommends locales && \
#     apt-get purge && \
#     apt-get clean && \
#     rm -rf /var/lib/apt/lists/*
ENV MAMBA_ROOT_PREFIX="/opt/conda"
WORKDIR "$MAMBA_ROOT_PREFIX"

# USER container_user
RUN chown -R ${NB_USER}:${NB_USER} ${HOME}
USER ${NB_USER}
SHELL ["conda", "run", "-n", "cenv", "/bin/bash", "-c"]

# RUN which pip
RUN pip install celltypist==1.3.0
############################ Setup: R packages #################################
USER micromamba
WORKDIR /tmp
ENTRYPOINT ["/bin/entrypoint.sh"]
CMD ["/bin/bash"]
################################################################################
