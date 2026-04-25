#!/usr/bin/env python
import os
import argparse
import scanpy as sc
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
    
    # Suppose adata is your AnnData object.
    # Define a function to check if a string starts with 'ENSG'
    def looks_like_ensembl(id_str):
        # A simple heuristic: Ensembl gene IDs typically start with 'ENSG'
        return isinstance(id_str, str) and id_str.startswith("ENSG")    

    # Get the gene symbols column as a Pandas Series
    gene_symbols_series = adata.var['gene_symbols'].astype(str)

    # Calculate the proportion of values that look like an Ensembl ID in the gene_symbols column.
    prop_in_gene_symbols = gene_symbols_series.str.startswith("ENSG").mean()

    # Calculate the proportion in the index (adata.var_names)
    prop_in_index = sum(looks_like_ensembl(gene) for gene in adata.var_names) / len(adata.var_names)

    print("Proportion of ENSG pattern in gene_symbols column: {:.2f}".format(prop_in_gene_symbols))
    print("Proportion of ENSG pattern in index: {:.2f}".format(prop_in_index))

    # Decide which one is which based on the proportions.
    # For example, if more than 80% of the values in a column start with 'ENSG', we assume that column contains Ensembl IDs.
    threshold = 0.8

    if prop_in_gene_symbols >= threshold and prop_in_index < threshold:
        # The gene_symbols column seems to be Ensembl IDs.
        # So we can assume the index holds the actual gene symbols.
        print("Using index as gene symbols and gene_symbols column as Ensembl IDs.")
        adata.var['gene_ids'] = gene_symbols_series
        # The index remains gene symbols (or is processed further if needed)
        
    elif prop_in_index >= threshold and prop_in_gene_symbols < threshold:
        # The index is mainly Ensembl IDs.
        # Therefore, the gene_symbols column should be used for gene symbols.
        print("Using gene_symbols column as gene symbols and index as Ensembl IDs.")
        adata.var['gene_ids'] = adata.var_names  # Preserve Ensembl IDs from the index
        # Now update the var_names to use the gene symbols
        adata.var_names = pd.Index(gene_symbols_series)
        
    elif prop_in_index >= threshold and prop_in_gene_symbols >= threshold:
        # Both columns have a high proportion of ENSG values.
        # This might indicate that both are Ensembl IDs or that there is an error.
        print("Warning: Both the index and gene_symbols appear to contain Ensembl IDs. A manual check may be needed.")
    else:
        # Neither column is clearly Ensembl IDs
        print("Warning: Neither the index nor gene_symbols contain a clear majority of Ensembl IDs. Please verify manually.")

    # Optionally, if the gene_symbols column is now redundant or incorrect, you might remove it:
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
