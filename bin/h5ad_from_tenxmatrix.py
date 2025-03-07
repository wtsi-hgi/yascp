#!/usr/bin/env python
import os
import argparse
import scanpy as sc
import anndata as ad
import pandas as pd
from scipy import io
import anndata as ad
import scipy.sparse as sp
    
# Define a no-op function (does nothing)
def no_op(*args, **kwargs):
    pass  # Do nothing

# Properly override AnnData class method
ad.AnnData.strings_to_categoricals = no_op  # This works for new instances


def mtx_to_h5ad(mtx_dir, output_file):
    # Paths to the input files
    adata = sc.read_10x_mtx(path=mtx_dir,
        var_names='gene_ids',
        make_unique=False
    )
    
    # Ensure feature metadata is Seurat-compatible
    if 'gene_symbols' not in adata.var.columns:
        adata.var['gene_symbols'] = adata.var.index.astype(str)

    if 'feature_types' not in adata.var.columns:
        adata.var['feature_types'] = 'Gene'

    if 'genome' not in adata.var.columns:
        adata.var['genome'] = 'GRCh38'  # Adjust if using different genome version

    # Ensure cell barcode metadata is explicitly stored
    if 'cell_barcode' not in adata.obs.columns:
        adata.obs['cell_barcode'] = adata.obs.index.astype(str)
    
    
    if any(adata.var.gene_symbols.str.contains('ENSG')):
        adata.var['gene_ids'] = adata.var.index
    else:
        adata.var_names =  pd.Index(adata.var.gene_symbols.astype(str))
        adata.var['gene_ids'] = adata.var['gene_symbols'] 
    
    del adata.var['gene_symbols']
    # Save to .h5ad file
    for col in adata.var.select_dtypes(include=['category']).columns:
        adata.var[col] = np.array(adata.var[col].astype(str), dtype=str)

    for col in adata.obs.select_dtypes(include=['category']).columns:
        adata.obs[col] = np.array(adata.obs[col].astype(str), dtype=str)
        


    if not sp.isspmatrix_csr(adata.X):
        adata.X = sp.csr_matrix(adata.X)
    
    adata.write(output_file)
    print(f"Data saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert Matrix Market files to .h5ad format")
    parser.add_argument("mtx_dir", help="Directory containing the matrix.mtx, features.tsv, and barcodes.tsv files")
    parser.add_argument("output_file", help="Output .h5ad file path")

    args = parser.parse_args()

    mtx_to_h5ad(args.mtx_dir, args.output_file)
