#!/usr/bin/env python3

__date__ = '2021-07-08'
__version__ = '0.0.1'

import argparse
import os
import pandas as pd
import glob
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
from distutils.version import LooseVersion
from typing import Dict
import tables
import scipy
import scipy.io
import scipy.sparse
import gzip
import pandas as pd
import numpy as np
import anndata
import scanpy as sc
import click
import logging

from plotnine import ggplot, geom_vline, geom_point, aes, labs, stat_smooth, facet_wrap, geom_histogram, geom_boxplot, theme, theme_bw, element_text, scale_x_continuous, scale_y_continuous, element_blank
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec

import tables
import numpy as np

def save_anndata_to_cellbender_h5(adata, path):
    with tables.open_file(path, mode='w') as f:
        # Transpose back to (n_barcodes, n_genes) as expected by CellBender
        X = adata.X.transpose().tocsc()
        f.create_array("/", 'data', X.data)
        f.create_array("/", 'indices', X.indices)
        f.create_array("/", 'indptr', X.indptr)
        f.create_array("/", 'shape', np.array(X.shape))  # should now be (n_cells, n_genes)

        # Barcodes
        f.create_array("/", 'barcodes', adata.obs_names.to_numpy(dtype='S'))

        # Gene metadata
        f.create_array("/", 'id', adata.var.index.to_numpy(dtype='S'))
        f.create_array("/", 'name', adata.var['gene_symbols'].to_numpy(dtype='S'))
        f.create_array("/", 'feature_type', adata.var['feature_type'].to_numpy(dtype='S'))

        # Write obs
        for key in adata.obs.columns:
            val = adata.obs[key].values
            if len(val) != adata.n_obs:
                print(f"Skipping obs key '{key}' due to length mismatch ({len(val)} vs {adata.n_obs})")
                continue
            try:
                f.create_array("/", key, np.asarray(val))
            except Exception as e:
                print(f"Could not write obs key '{key}': {e}")

        # Write obsm
        for key in adata.obsm_keys():
            val = adata.obsm[key]
            if val.shape[0] != adata.n_obs:
                print(f"Skipping obsm key '{key}' due to mismatch ({val.shape[0]} vs {adata.n_obs})")
                continue
            try:
                f.create_array("/", key, val)
            except Exception as e:
                print(f"Could not write obsm key '{key}': {e}")

        # Write uns
        for key, val in adata.uns.items():
            try:
                f.create_array("/", key, np.atleast_1d(val))
            except Exception as e:
                print(f"Could not write uns key '{key}': {e}")


def dict_from_h5(file: str) -> Dict[str, np.ndarray]:
    """Read in everything from an h5 file and put into a dictionary."""
    d = {}
    with tables.open_file(file) as f:
        # read in everything
        for array in f.walk_nodes("/", "Array"):
            d[array.name] = array.read()
    return d


def anndata_from_h5(
    file: str,
    analyzed_barcodes_only: bool = True
 ) -> 'anndata.AnnData':
    """Load an output h5 file into an AnnData object for downstream work.

    Args:
        file: The h5 file
        analyzed_barcodes_only: False to load all barcodes, so that the size of
            the AnnData object will match the size of the input raw count
            matrix. True to load a limited set of barcodes: only those
            analyzed by the algorithm. This allows relevant latent
            variables to be loaded properly into adata.obs and adata.obsm,
            rather than adata.uns.

    Returns:
        adata: The anndata object, populated with inferred latent variables
            and metadata.
    """
    d = dict_from_h5(file)
    X = scipy.sparse.csc_matrix(
        (d.pop('data'), d.pop('indices'), d.pop('indptr')),
        shape=d.pop('shape')
    ).transpose().tocsr()

    if analyzed_barcodes_only:
        if 'barcodes_analyzed_inds' in d.keys():
            X = X[d['barcodes_analyzed_inds'], :]
            d['barcodes'] = d['barcodes'][d['barcodes_analyzed_inds']]
        elif 'barcode_indices_for_latents' in d.keys():
            X = X[d['barcode_indices_for_latents'], :]
            d['barcodes'] = d['barcodes'][d['barcode_indices_for_latents']]
        else:
            print(
                'Warning: analyzed_barcodes_only=True, but the key ',
                '"barcodes_analyzed_inds" or "barcode_indices_for_latents" ',
                'is missing from the h5 file. ',
                'Will output all barcodes, and proceed as if ',
                'analyzed_barcodes_only=False'
            )

    print(d.keys())

    # Construct the count matrix.
    if 'gene_names' in d.keys():
        gene_symbols = d.pop('gene_names').astype(str)
    else:
        gene_symbols = d.pop('name').astype(str)
    adata = anndata.AnnData(
        X=X,
        obs={'barcode': d.pop('barcodes').astype(str)},
        var={
            'gene_ids': d.pop('id').astype(str),
            'gene_symbols': gene_symbols,
            'feature_type':d.pop('feature_type').astype(str),
        }
    )
    adata = adata[:,adata.var.query('feature_type=="Gene Expression"').index]
    adata.obs.set_index('barcode', inplace=True)
    adata.var.set_index('gene_ids', inplace=True)
    
    
    
    # Add other information to the adata object in the appropriate slot.
    for key, value in d.items():
        try:
            value = np.asarray(value)
            if len(value.shape) == 0:
                adata.uns[key] = value
            elif value.shape[0] == X.shape[0]:
                if (len(value.shape) < 2) or (value.shape[1] < 2):
                    adata.obs[key] = value
                else:
                    adata.obsm[key] = value
            elif value.shape[0] == X.shape[1]:
                if value.dtype.name.startswith('bytes'):
                    adata.var[key] = value.astype(str)
                else:
                    adata.var[key] = value
            else:
                adata.uns[key] = value
        except Exception:
            print(
                'Unable to load data into AnnData: ', key, value, type(value)
            )

    if analyzed_barcodes_only:
        cols = adata.obs.columns[
            adata.obs.columns.str.startswith('barcodes_analyzed')
            | adata.obs.columns.str.startswith('barcode_indices')
        ]
        for col in cols:
            try:
                del adata.obs[col]
            except Exception:
                pass

    if adata.obs_names.duplicated().any():
        print("Warning: duplicated barcodes found — making obs_names unique.")
        adata.obs_names_make_unique()
    return adata


def remap_barcodes(adata, mapping_tsv_path):
    """
    Replace barcodes in adata.obs_names and adata.obs['barcode'] using a mapping file.
    
    Parameters:
        adata: AnnData object
        mapping_tsv_path: path to TSV file with columns 'original_barcode' and 'encoded_barcode'
    """
    # Load the mapping
    mapping_df = pd.read_csv(mapping_tsv_path, sep='\t')
    mapping = dict(zip(mapping_df['encoded_barcode'], mapping_df['original_barcode']))

    # Check how many barcodes are matched
    matched = adata.obs_names.isin(mapping.keys()).sum()
    print(f"Replacing {matched} barcodes using provided mapping.")

    # Remap
    new_barcodes = adata.obs_names.map(mapping)
    # If any aren't mapped, warn
    if new_barcodes.isna().any():
        unmapped = adata.obs_names[new_barcodes.isna()]
        print(f"Warning: {unmapped.shape[0]} barcodes not found in mapping.")
        # Optional: keep original names for those or raise an error
        new_barcodes = new_barcodes.fillna(adata.obs_names)
    # Assign new barcodes
    adata.obs_names = new_barcodes
    return adata

mapping_file = "./txd_input/barcode_mapping.tsv"
if os.path.exists(mapping_file):
    all_files_to_resave = glob.glob("./*.h5")
    for f1 in all_files_to_resave:
        adata_cellbender = anndata_from_h5(f1,analyzed_barcodes_only=False)
        adata_cellbender2 = remap_barcodes(
            adata_cellbender,
            "./txd_input/barcode_mapping.tsv"
        )
        save_anndata_to_cellbender_h5(adata_cellbender, f1)

print('Done')