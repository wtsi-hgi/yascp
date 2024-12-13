#!/usr/bin/env python

__date__ = '2020-05-26'
__version__ = '0.0.1'

import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec
import numpy as np
import pandas as pd
import scipy as sci
import scanpy as sc
import copy
import scvi
import argparse
import glob
import pandas as pd



# This code will gather the most important graphs and place them in a summary folder.
parser = argparse.ArgumentParser(
    description="""
        Collect inputs
        """
)

parser.add_argument(
        '-h5ad_file', '--h5ad_file',
        action='store',
        dest='h5ad_file',
        required=True,
        help='Input h5ad_file after normalisation'
)



options = parser.parse_args()

sc.set_figure_params(figsize=(5,5), dpi=150)

#import single cell data and CITE-seq data
# SLEmap = sc.read('adata-normalized.h5ad')
SLEmap = sc.read(options.h5ad_file, backed ='r')

all_cite_files = glob.glob("./*/*.matrix.csv")
CITE = pd.DataFrame()
for f1 in all_cite_files:
    c1 = pd.read_csv(f1,index_col=0)
    CITE=pd.concat([CITE,c1])

# merge citeseq data to obs of single cell data 
SLEmap.obs["Barcode"] = SLEmap.obs.index.str.split('__').str[0]
SLEmap.obs = SLEmap.obs.merge(CITE, left_on=['Barcode'],right_index=True,how='left', indicator=True)


# make citeseq data to obsm single cell data : required for totalVI
CITE_2 = SLEmap.obs[CITE.columns].copy()
SLEmap.obsm['protein_expression'] = CITE_2

# keep only cells passing QC and highly variable genes
SLEmap = SLEmap[
    (SLEmap.obs["cell_passes_qc"] & SLEmap.obs["cell_passes_hard_filters"]),
    SLEmap.var["highly_variable"]
]

#run totalVI
SLEmap = SLEmap.to_memory().copy()
scvi.model.TOTALVI.setup_anndata(SLEmap, protein_expression_obsm_key="protein_expression",batch_key="experiment_id")
model = scvi.model.TOTALVI(SLEmap,latent_distribution="normal",n_layers_decoder=2)

model.train()
model.save("./scvi_model",adata=SLEmap, overwrite=True)


# Get latent expression from model: used for UMAP calculations
SLEmap.obsm["X_totalVI"] = model.get_latent_representation()


# Get denoised protein values: can fail with low memory: use 400 minimum for n_samples=25. 
# If it keeps failing lower n_samples or don't use it and use dsb values instead
rna, protein = model.get_normalized_expression(n_samples=25,return_mean=True)
SLEmap.obsm["denoised_protein"] = protein

sc.pp.neighbors(SLEmap, use_rep="X_totalVI",n_neighbors=15)
sc.tl.umap(SLEmap)


sc.set_figure_params(figsize=(5,5), dpi=150)

def panel_grid(hspace, wspace, ncols, num_panels):
    """Init plot."""
    n_panels_x = min(ncols, num_panels)
    n_panels_y = np.ceil(num_panels / n_panels_x).astype(int)
    if wspace is None:
        #  try to set a wspace that is not too large or too small given the
        #  current figure size
        wspace = 0.75 / rcParams['figure.figsize'][0] + 0.02
    # each panel will have the size of rcParams['figure.figsize']
    fig = plt.figure(
        figsize=(
            n_panels_x * rcParams['figure.figsize'][0] * (1 + wspace),
            n_panels_y * rcParams['figure.figsize'][1],
        )
    )
    left = 0.2 / n_panels_x
    bottom = 0.13 / n_panels_y
    gs = gridspec.GridSpec(
        nrows=n_panels_y,
        ncols=n_panels_x,
        left=left,
        right=1 - (n_panels_x - 1) * left - 0.01 / n_panels_x,
        bottom=bottom,
        top=1 - (n_panels_y - 1) * bottom - 0.1 / n_panels_y,
        hspace=hspace,
        wspace=wspace
    )
    return fig, gs

fig, grid = panel_grid(
    hspace=0.125*2,
    wspace=None,
    ncols=4,
    num_panels=2
)

sc.pl.umap(
    SLEmap,
    color=["Azimuth:predicted.celltype.l1","Azimuth:predicted.celltype.l2"],
    # Setting a smaller point size to get prevent overlap
    size=0.7,save='plots.pdf'
)

fig.savefig(
    'umap.png',
    #dpi=300,
    bbox_inches='tight'
)

SLEmap.write("./totalVI_integrated.h5ad")