#!/usr/bin/env python


__date__ = '2021-01-20'
__version__ = '0.0.1'
# sklearn version used = 0.24.2
import argparse
from distutils.version import LooseVersion
import numpy as np
import scanpy as sc
import pandas as pd
import glob
       
def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Performs automatic outlier detection over an anndata dataset,
            plotting the outlier cells.
            """
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file.'
    )

    parser.add_argument(
        '--anndata_compression_opts',
        action='store',
        dest='anndata_compression_opts',
        default=4,
        type=int,
        help='Compression level in anndata. A larger value decreases disk \
            space requirements at the cost of compression time. \
            (default: %(default)s)'
    )

    options = parser.parse_args()

    # Load the AnnData file.
    adata = sc.read_h5ad(filename=options.h5)
    if 'cell_passes_qc' in adata.obs.columns:
        del adata.obs['cell_passes_qc']
        
    files = glob.glob('./outlier_filtered_adata-outliers_filtered__*.tsv')
    combo_files = pd.DataFrame()
    for f1 in files:
        f2= pd.read_csv(f1,sep='\t',index_col=0)
        combo_files = pd.concat([combo_files,f2],axis=1)
        
    adata.obs = pd.concat([adata.obs,combo_files],axis=1)

    adata.write(
        '4.outlier_filtered_adata.h5ad',
        compression='gzip',
        compression_opts=options.anndata_compression_opts
    )


if __name__ == '__main__':
    main()
