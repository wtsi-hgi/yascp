#!/usr/bin/env python3

__date__ = '2021-09-29'
__version__ = '0.0.1'
# for help, run: python3 run_celltypist.py --help
# command with example inputs arguments:
# python3 run_celltypist.py.py \
#  --samplename Pilot_study_of_dissociation_methods_for_human_gut_tissues8024876 \
#  --filtered_matrix_h5 filtered_feature_bc_matrix.h5 \
#  --celltypist_model Immune_Blood_Low.pkl \
#  --output_dir $PWD/outputs

# https://pypi.org/project/celltypist/

# import python libraries:
import logging
import click
import sys 
import argparse
import os
import csv
import random
import numpy as np
import pandas as pd
import pandas
import os
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
os.environ['HOME'] = '/tmp'
import scanpy as sc
import celltypist
from celltypist import models
from pathlib import Path
import shutil

total, used, free = shutil.disk_usage("/tmp")

logging.info("Total: %d GiB" % (total // (2**30)))
logging.info("Used: %d GiB" % (used // (2**30)))
logging.info("Free: %d GiB" % (free // (2**30)))
# CLI arguments:
@click.command()

# required arguments:
@click.option('--samplename', required=True, type=str,
              help='sample name of cellranger experiment deconvoluted by Vireo')

@click.option('--filtered_matrix_h5', required=True, type=click.Path(exists=True),
              help='path to cellranger output filtered h5, which probably ends with .filtered_feature_bc_matrix.h5')

@click.option('-m','--celltypist_model', required=True, type=str,
              help='celltypist model .pkl from models.download_models(force_update = True)')

@click.option('-o','--output_dir', required=True, type=str,
              help='output directory for script output plots and files')

@click.option('-c','--anndata_compression_level', default=6,
              type=click.IntRange(1, 9, clamp=True), show_default=True,
              help='Gzip compression level for scanpy write of AnnData hdf5 objects. Integer in range 1 to 9')

@click.option('-g','--input_h5_genome_version', default="GRCh38", show_default=True, type=str,
              help='True or False: whether to write donor level scanpy hdf5 objects to dir --output_dir')

# Optional arguments:
@click.option('-p','--sample_plot_probs', is_flag=True, default=False, type=bool,
              help='True or False: whether or not to plot probabilities per cell types and per sample')

def run_celltypist(samplename, filtered_matrix_h5, celltypist_model,
                   output_dir, anndata_compression_level,input_h5_genome_version, sample_plot_probs):
    """process cellranger output filtered h5 so that it can be fed to Celltypist"""
    logging.info('running run_celltypist() function..')

    # Set seed for reproducibility
    seed_value = 0
    # 0. Set `PYTHONHASHSEED` environment variable at a fixed value
    # os.environ['PYTHONHASHSEED']=str(seed_value)
    # 1. Set `python` built-in pseudo-random generator at a fixed value
    random.seed(seed_value)
    # 2. Set `numpy` pseudo-random generator at a fixed value
    np.random.seed(seed_value)
    fixed_h5 = 'fixed_genome.h5'
    try:
        adata = sc.read_10x_h5(filtered_matrix_h5)
    except:
        try:
            # We are loading h5ad instad of h5
            adata = sc.read_h5ad(filtered_matrix_h5)
        except:
            # h5 file may have amissing genome version
            import tables
            logging.info('fixing orig_h5')
            orig_h5 = filtered_matrix_h5
            fixed_h5 = 'fixed_genome.h5'
            tables.copy_file(orig_h5, fixed_h5, overwrite = True)
            with tables.open_file(fixed_h5, "r+") as f:
                n = f.get_node("/matrix/features")
                n_genes = f.get_node("/matrix/shape")[0]
                if "genome" not in n:
                    f.create_array(n, "genome", np.repeat(input_h5_genome_version, n_genes))
            # read-in cellranger 10x data produced by 'cellranger count':
            logging.info('fixed orig_h5 into fixed_h5')
            adata = sc.read_10x_h5(fixed_h5) #, genome='background_removed')
        
    try:
        os.remove(fixed_h5)
        logging.info('removed fixed h5 file for memory saving - should really fix it in the cellbender pipeline')
    except:
        logging.info('no file to remove')
    logging.info('loadin sc.read_10x_h5() done.')
    logging.info(adata.var)
    # adata.var['ENSG']=adata.var.index
    # adata.var.index=adata.var['gene_symbols']
    logging.info(adata.obsm)
    logging.info("n cells len(adata.obs): " + str(len(adata.obs)))
    # logging.info('exiting..'); sys.exit()

    # how to prep h5 data for celltypist,  
    # from README https://github.com/Teichlab/celltypist : 
    # CellTypist requires a logarithmised and normalised expression matrix stored in the AnnData
    #   (log1p normalised expression to 10,000 counts per cell).
    #   CellTypist will try the .X attribute first, and if it does not suffice,
    #       try the .raw.X attribute.
    # If none of them fit into the desired data type or the expression matrix is not properly normalised, an error will be raised.
    logging.info('... running sc.pp.normalize_total(adata, target_sum=1e4)')
    sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
    logging.info('... running sc.pp.log1p(adata, copy = False)')
    sc.pp.log1p(adata, copy = False)

    # output_file = output_dir + '/' + output_filename_prefix + filtered_matrix_h5
    # logging.info('output file: ' + output_file)
    # adata.write_h5ad(filename=output_file, compression='gzip',
    #                  compression_opts= anndata_compression_level)

    # check The expression matrix (adata_2000.X) is pre-processed
    #  (and required) as log1p normalised expression to 10,000 counts per cell
        
    # https://www.celltypist.org/tutorials
    # https://colab.research.google.com/github/Teichlab/celltypist/blob/main/notebook/celltypist_tutorial.ipynb#scrollTo=ultimate-pilot
    # Enabling `force_update = True` will overwrite existing (old) models.

    logging.info("celltypist_model: " + celltypist_model)
    celltypist_model1 = os.path.splitext(os.path.basename(celltypist_model))[0]
    if not os.path.isfile(celltypist_model):
        p1 =  '/'.join(__file__.split('/')[:-1])
        celltypist_model=f"{p1}/../assets/celltypist/{celltypist_model.split('/')[-1]}"
    model = models.Model.load(model = celltypist_model ) # model = 'Immune_All_Low.pkl')
    model.description

    # Not run; predict cell identities using this loaded model.
    # predictions = celltypist.annotate(adata_2000, model = model, majority_voting = True)   
    # Alternatively, just specify the model name (recommended as this ensures the model is intact every time it is loaded).
    logging.info("... running celltypist.annotate(adata, model = model, majority_voting = True)")
    try:
        predictions = celltypist.annotate(adata, model = model, majority_voting = True, transpose_input = True)
    except:
        adata.var['gene_ids'] = adata.var.index
        adata.var_names =  pd.Index(adata.var.gene_symbols.astype(str))
        predictions = celltypist.annotate(adata, model = model, majority_voting = True, transpose_input = True)
    
    # By default (majority_voting = False), CellTypist will infer the identity of
    # each query cell independently. This leads to raw predicted cell type labels,
    # and usually finishes within seconds or minutes depending on the size of the query data.
    # You can also turn on the majority-voting classifier (majority_voting = True),
    # which refines cell identities within local subclusters after an over-clustering approach
    # at the cost of increased runtime.

    logging.info("... predictions.predicted_labels:")
    logging.info(predictions.predicted_labels)

    # Export the three results to csv tables.
    logging.info("... predictions.to_table in folder " + output_dir)
    
    adata = predictions.to_adata()

    predictions.to_table(folder = output_dir, prefix = samplename + '___'+celltypist_model1+'___')
    Data = adata.obs
    try:
        Data= Data.drop('cell_barcode',axis=1)
    except:
        _=''
        
    Data=Data.add_prefix(f'{celltypist_model1}:')
    Data.to_csv(f'{output_dir}/{samplename}___{celltypist_model1}___predicted_labels.csv')
    ###predictions.to_table(folder = os.getcwd())
    logging.info("... predictions.to_plots")
    # Visualise the predicted cell types overlaid onto the UMAP.
    predictions.to_plots(folder = output_dir, prefix = samplename + '_')
    ###predictions.to_plots(folder = os.getcwd())
    # Visualise the decision scores and probabilities of each cell type overlaid onto the UMAP as well.
    if sample_plot_probs:
        folder_plot_probs = output_dir + '/plot_prob'
        if not os.path.exists(folder_plot_probs):
            os.makedirs(folder_plot_probs)
        predictions.to_plots(folder = folder_plot_probs, prefix = samplename + '_prob_', plot_probability = True)

    logging.info("... script run_celltypist() done.")
    
if __name__ == '__main__':
    # set logging level and handler:
    logging.basicConfig(level=logging.INFO,
                        format="%(asctime)s [%(levelname)s] %(message)s",
                        handlers=[logging.StreamHandler()]) # logging.FileHandler("debug.log"),
    run_celltypist()

