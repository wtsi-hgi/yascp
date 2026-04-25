#!/usr/bin/env python
import os
import argparse
import scanpy as sc
import anndata as ad
import pandas as pd
from scipy import io

def mtx_to_h5ad(mtx_dir, output_file):
    # Paths to the input files
    adata = sc.read_10x_mtx(path=mtx_dir,
        var_names='gene_ids',
        make_unique=False
    )
    
    # Save to .h5ad file
    adata.write(output_file)
    print(f"Data saved to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert Matrix Market files to .h5ad format")
    parser.add_argument("mtx_dir", help="Directory containing the matrix.mtx, features.tsv, and barcodes.tsv files")
    parser.add_argument("output_file", help="Output .h5ad file path")

    args = parser.parse_args()

    mtx_to_h5ad(args.mtx_dir, args.output_file)
