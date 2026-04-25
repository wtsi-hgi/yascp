#!/usr/bin/env python
import anndata
import scipy.io
import os
import pandas as pd
import numpy as np
import gzip
from distutils.version import LooseVersion
import argparse

def h5ad_to_tenxmatrix(
    adata,
    out_file='',
    out_dir='tenx_from_adata',
    verbose=True
):
    """Convert AnnData object to 10x-like Matrix Market format.

    Parameters
    ----------
    adata : anndata.AnnData
        Annotated data matrix.
    out_file : string
        Prefix for output files.
    out_dir : string
        Directory for output files.
    verbose : boolean
        Verbosity of the function.

    Returns
    -------
    execution_code : int
    """
    # Make the output directory if it does not exist.
    if out_dir == '':
        out_dir = os.getcwd()
    else:
        os.makedirs(out_dir, exist_ok=True)

    # Set up out_file
    if out_file != '':
        out_file = '{}-'.format(out_file)

    # Get compression opts for pandas
    compression_opts = 'gzip'
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(method='gzip', compresslevel=9)

    # First filter out any cells that have 0 total counts
    zero_count_cells = adata.obs_names[np.where(adata.X.sum(axis=1) == 0)[0]]
    if verbose:
        print("Filtering {}/{} cells with 0 counts.".format(
            len(zero_count_cells),
            adata.n_obs
        ))
    adata = adata[adata.obs_names.difference(zero_count_cells, sort=False)]

    # Save the barcodes.
    out_f = os.path.join(
        out_dir,
        '{}barcodes.tsv.gz'.format(out_file)
    )
    if verbose:
        print('Writing {}'.format(out_f))
    pd.DataFrame(adata.obs.index).to_csv(
        out_f,
        sep='\t',
        compression=compression_opts,
        index=False,
        header=False
    )

    # Save the features.
    out_f = os.path.join(
        out_dir,
        '{}features.tsv.gz'.format(out_file)
    )
    if verbose:
        print('Writing {}'.format(out_f))

    # Check if read in by gene_symbol or gene_id
    if 'gene_ids' in adata.var.columns:
        gene_var = 'gene_ids'
    elif 'gene_symbols' in adata.var.columns:
        gene_var = 'gene_symbols'
    else:
        raise Exception(
            'Could not find "gene_symbols" or "gene_ids" in adata.var'
        )
    df_features = pd.DataFrame(
        data=[
            adata.var.index.values,
            adata.var.loc[:, gene_var].values,
            ["Gene Expression"] * adata.n_vars
        ]
    ).T

    df_features.to_csv(
        out_f,
        sep='\t',
        compression=compression_opts,
        index=False,
        header=False
    )

    # Save the count-adjusted matrix
    out_mtx = adata.X.transpose()
    if not isinstance(out_mtx, scipy.sparse.csr.csr_matrix):
        out_mtx = scipy.sparse.csr_matrix(out_mtx)
    out_f = os.path.join(
        out_dir,
        '{}matrix.mtx.gz'.format(out_file)
    )
    if verbose:
        print('Writing {}'.format(out_f))
    with gzip.open(out_f, 'wb', compresslevel=9) as f:
        scipy.io.mmwrite(
            f,
            out_mtx,
            comment='metadata_json: {"format_version": 2}'
        )

    if verbose:
        print('Done.')

    return 0

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert .h5ad file to 10x-like Matrix Market format")
    parser.add_argument("h5ad_path", help="Path to the .h5ad file")
    parser.add_argument("output_dir", help="Directory to save the output files")
    parser.add_argument("--out_file", default="", help="Prefix for output files")
    parser.add_argument("--verbose", action="store_true", help="Increase output verbosity")

    args = parser.parse_args()

    adata = anndata.read_h5ad(args.h5ad_path)
    h5ad_to_tenxmatrix(adata, args.out_file, args.output_dir, args.verbose)
