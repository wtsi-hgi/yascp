#!/usr/bin/env python

__date__ = '2020-04-03'
__version__ = '0.0.1'
import numpy as np
import pandas as pd
import scanpy as sc
import argparse
import os
import yaml
from distutils.version import LooseVersion
import csv
y_n_print=False
def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Filter and merge h5ad data. Save to AnnData object.
            """
    )
    
    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )    
    
    parser.add_argument(
        '-ad', '--h5addata_file',
        action='store',
        dest='adata_file',
        required=True,
        help='adata_file'
    )
 
    parser.add_argument(
        '-ncpu', '--number_cpu',
        action='store',
        dest='ncpu',
        default=2,
        type=int,
        help='Number of CPUs to use.\
            (default: %(default)s)'
    ) 
    
    parser.add_argument(
        '-name', '--experiment_name',
        action='store',
        dest='experiment_name',
        default='all_samples',
        help='experiment name to determine whether this sample has to be treated differently'
    )

    parser.add_argument(
        '-pyml', '--params_yaml',
        action='store',
        dest='pyml',
        default='',
        help='YAML file containing cell filtering and downsampling (e.g.,\
            of total number of cells or reads per cell) parameters.\
            If file is not provided, no filtering or downsampling is\
            performed.\
            (default: %(default)s)'
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
    
    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='adata',
        help='Basename of output anndata file, assuming output in current \
            working directory. Will have .h5ad appended.\
            (default: %(default)s)'
    )
        
    parser.add_argument(
        '-mk', '--metadata_key',
        action='store',
        dest='mk',
        default='experiment_id',
        help='Key to link metadata to h5addata_file experiment_id column.\
            (default: %(default)s)'
    )    
    
    options = parser.parse_args()
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    sc.settings.n_jobs = options.ncpu  # number CPUs

    # Read in the parameters for downsampling and cell filters.
    if options.pyml == '':
        params_dict = {}
    else:
        with open(options.pyml, 'r') as f:
            params_dict = yaml.safe_load(f)
        params_dict = params_dict['sample_qc']
 
    anndata_compression_opts=options.anndata_compression_opts
    compression_opts = 'gzip'
    
    if LooseVersion(pd.__version__) > '1.0.0':
        compression_opts = dict(
            method='gzip',
            compresslevel=anndata_compression_opts
        )         

        
    params_filters_check = [
        'cell_filters',
        'downsample_cells_fraction',
        'downsample_cells_n',
        'downsample_feature_counts'
    ]
    
    for i in params_filters_check:
        if i not in params_dict:
            if i != 'cell_filters':
                params_dict[i] = {'value': ''}
    # Check for validity of filters specified in params_dict
    param_filters_check = [
        'downsample_cells_fraction',
        'downsample_cells_n'
    ]
    if all(params_dict[k]['value'] != '' for k in param_filters_check):
        raise Exception(
            'Error check the params. Both {} and {} are set.'.format(
                'downsample_cells_fraction',
                'downsample_cells_n'
            )
        )       
    ad_file = options.adata_file
    adata = sc.read_h5ad(ad_file)

    adata_orig_cols = list(adata.obs.columns)
    adata_orig_cols.append("donor")

    # Record the total number of cells for this experiment_id
    n_cells_dict = {}
    n_cells_dict[options.experiment_name] = {}
    n_cells_dict[options.experiment_name]['before_filters'] = adata.n_obs

    # Label mitochondrial encoded trancripts
    # This includes:
    # * Mitochondrially encoded protein coding genes
    # * Mitochondrially encoded transcribed regions
    # * Mitochondrially encoded RNAs
    #
    # The below gene list was downloaded from on 3 Aug 2020:
    # https://www.genenames.org/data/genegroup/#!/group/1972
    gene_group__mito_transcript = [
        'MT-7SDNA',
        'MT-ATP6',
        'MT-ATP8',
        'MT-ATT',
        'MT-CO1',
        'MT-CO2',
        'MT-CO3',
        'MT-CSB1',
        'MT-CSB2',
        'MT-CSB3',
        'MT-CYB',
        'MT-HPR',
        'MT-HSP1',
        'MT-HSP2',
        'MT-LIPCAR',
        'MT-LSP',
        'MT-ND1',
        'MT-ND2',
        'MT-ND3',
        'MT-ND4',
        'MT-ND4L',
        'MT-ND5',
        'MT-ND6',
        'MT-OHR',
        'MT-OLR',
        'MT-RNR1',
        'MT-RNR2',
        'MT-RNR3',
        'MT-TA',
        'MT-TAS',
        'MT-TC',
        'MT-TD',
        'MT-TE',
        'MT-TER',
        'MT-TF',
        'MT-TFH',
        'MT-TFL',
        'MT-TFX',
        'MT-TFY',
        'MT-TG',
        'MT-TH',
        'MT-TI',
        'MT-TK',
        'MT-TL1',
        'MT-TL2',
        'MT-TM',
        'MT-TN',
        'MT-TP',
        'MT-TQ',
        'MT-TR',
        'MT-TS1',
        'MT-TS2',
        'MT-TT',
        'MT-TV',
        'MT-TW',
        'MT-TY'
    ]
    adata.var['gene_group__mito_transcript'] = [
        x in gene_group__mito_transcript for x in adata.var['gene_symbols']
    ]
    # mito_gene_list = sc.queries.mitochondrial_genes() # another query
    # adata.var['gene_group__mito_transcript'] = [
    #     x.startswith('MT-') for x in adata.var['gene_symbols']
    # ]
    # use this if var_names='gene_symbols' in sc.read_h5ad_mtx
    # adata.var['mito_gene'] = [
    #     x.startswith('MT-') for x in adata.var_names
    # ]

    # Label mitochondrial encoded proteins
    #
    # The below gene list was downloaded from on 3 Aug 2020:
    # https://www.genenames.org/data/genegroup/#!/group/1974
    gene_group__mito_protein = [
        'MT-ATP6',
        'MT-ATP8',
        'MT-CO1',
        'MT-CO2',
        'MT-CO3',
        'MT-CYB',
        'MT-ND1',
        'MT-ND2',
        'MT-ND3',
        'MT-ND4',
        'MT-ND4L',
        'MT-ND5',
        'MT-ND6'
    ]
    adata.var['gene_group__mito_protein'] = [
        x in gene_group__mito_protein for x in adata.var['gene_symbols']
    ]

    # Label ribosomal protein genes
    # Ribosomal protein: A ribosomal protein is any of
    # the proteins that, in conjunction with rRNA, make up the ribosomal
    # subunits involved in the cellular process of translation. A large
    # part of the knowledge about these organic molecules has come from the
    # study of E. coli ribosomes. Most ribosomal proteins have been
    # isolated and specific antibodies have been produced. These, together
    # with electronic microscopy and the use of certain reactives, have
    # allowed for the determination of the topography of the proteins in
    # ribosome. E.coli, other bacteria and Archaea have a 30S small subunit
    # and a 50S large subunit, whereas humans and yeasts have a 40S small
    # subunit and a 60S large subunit. Equivalent subunits are frequently
    # numbered differently between bacteria, Archaea, yeasts and humans.
    #
    # The below gene list was downloaded from on 3 Aug 2020:
    # https://www.genenames.org/data/genegroup/#!/group/1054
    gene_group__ribo_protein = set([
        'DAP3',
        'FAU',
        'MRPL1',
        'MRPL1',
        'MRPL10',
        'MRPL10',
        'MRPL11',
        'MRPL11',
        'MRPL12',
        'MRPL12',
        'MRPL13',
        'MRPL13',
        'MRPL14',
        'MRPL14',
        'MRPL15',
        'MRPL15',
        'MRPL16',
        'MRPL16',
        'MRPL17',
        'MRPL17',
        'MRPL18',
        'MRPL18',
        'MRPL19',
        'MRPL19',
        'MRPL2',
        'MRPL2',
        'MRPL20',
        'MRPL20',
        'MRPL21',
        'MRPL21',
        'MRPL22',
        'MRPL22',
        'MRPL23',
        'MRPL23',
        'MRPL24',
        'MRPL24',
        'MRPL27',
        'MRPL27',
        'MRPL28',
        'MRPL28',
        'MRPL3',
        'MRPL3',
        'MRPL30',
        'MRPL30',
        'MRPL32',
        'MRPL32',
        'MRPL33',
        'MRPL33',
        'MRPL34',
        'MRPL34',
        'MRPL35',
        'MRPL35',
        'MRPL36',
        'MRPL36',
        'MRPL37',
        'MRPL37',
        'MRPL38',
        'MRPL38',
        'MRPL39',
        'MRPL39',
        'MRPL4',
        'MRPL4',
        'MRPL40',
        'MRPL40',
        'MRPL41',
        'MRPL41',
        'MRPL42',
        'MRPL42',
        'MRPL43',
        'MRPL43',
        'MRPL44',
        'MRPL44',
        'MRPL45',
        'MRPL45',
        'MRPL46',
        'MRPL46',
        'MRPL47',
        'MRPL47',
        'MRPL48',
        'MRPL48',
        'MRPL49',
        'MRPL49',
        'MRPL50',
        'MRPL50',
        'MRPL51',
        'MRPL51',
        'MRPL52',
        'MRPL52',
        'MRPL53',
        'MRPL53',
        'MRPL54',
        'MRPL54',
        'MRPL55',
        'MRPL55',
        'MRPL57',
        'MRPL57',
        'MRPL58',
        'MRPL9',
        'MRPS10',
        'MRPS10',
        'MRPS11',
        'MRPS11',
        'MRPS12',
        'MRPS12',
        'MRPS14',
        'MRPS14',
        'MRPS15',
        'MRPS15',
        'MRPS16',
        'MRPS16',
        'MRPS17',
        'MRPS17',
        'MRPS18A',
        'MRPS18A',
        'MRPS18B',
        'MRPS18B',
        'MRPS18C',
        'MRPS18C',
        'MRPS2',
        'MRPS2',
        'MRPS21',
        'MRPS21',
        'MRPS22',
        'MRPS22',
        'MRPS23',
        'MRPS23',
        'MRPS24',
        'MRPS24',
        'MRPS25',
        'MRPS25',
        'MRPS26',
        'MRPS26',
        'MRPS27',
        'MRPS27',
        'MRPS28',
        'MRPS28',
        'MRPS30',
        'MRPS30',
        'MRPS31',
        'MRPS31',
        'MRPS33',
        'MRPS33',
        'MRPS34',
        'MRPS34',
        'MRPS35',
        'MRPS35',
        'MRPS36',
        'MRPS36',
        'MRPS5',
        'MRPS6',
        'MRPS7',
        'MRPS9',
        'RPL10',
        'RPL10A',
        'RPL10L',
        'RPL11',
        'RPL12',
        'RPL13A',
        'RPL14',
        'RPL15',
        'RPL17',
        'RPL18A',
        'RPL19',
        'RPL21',
        'RPL22',
        'RPL23',
        'RPL23A',
        'RPL24',
        'RPL26',
        'RPL26L1',
        'RPL27',
        'RPL27A',
        'RPL28',
        'RPL29',
        'RPL3',
        'RPL30',
        'RPL31',
        'RPL32',
        'RPL34',
        'RPL35',
        'RPL35A',
        'RPL36',
        'RPL36A',
        'RPL36AL',
        'RPL37',
        'RPL37A',
        'RPL38',
        'RPL39',
        'RPL39L',
        'RPL3L',
        'RPL4',
        'RPL41',
        'RPL5',
        'RPL6',
        'RPL7',
        'RPL7A',
        'RPL7L1',
        'RPL8',
        'RPL9',
        'RPLP0',
        'RPLP1',
        'RPLP2',
        'RPS10',
        'RPS11',
        'RPS12',
        'RPS13',
        'RPS14',
        'RPS15',
        'RPS15A',
        'RPS16',
        'RPS17',
        'RPS18',
        'RPS19',
        'RPS2',
        'RPS20',
        'RPS21',
        'RPS23',
        'RPS24',
        'RPS25',
        'RPS26',
        'RPS27',
        'RPS27A',
        'RPS27L',
        'RPS28',
        'RPS29',
        'RPS3',
        'RPS3A',
        'RPS4X',
        'RPS4Y1',
        'RPS4Y2',
        'RPS5',
        'RPS6',
        'RPS7',
        'RPS8',
        'RPS9',
        'UBA52'
    ])
    adata.var['gene_group__ribo_protein'] = [
        x in gene_group__ribo_protein for x in adata.var['gene_symbols']
    ]

    # Label ribosomal RNA
    #
    # The below gene list was downloaded from on 3 Aug 2020:
    # https://www.genenames.org/data/genegroup/#!/group/848
    gene_group__ribo_rna = [
        'MT-RNR1',
        'MT-RNR2',
        'RNA18S1',
        'RNA18S2',
        'RNA18S3',
        'RNA18S4',
        'RNA18S5',
        'RNA18SN1',
        'RNA18SN2',
        'RNA18SN3',
        'RNA18SN4',
        'RNA18SN5',
        'RNA28S1',
        'RNA28S2',
        'RNA28S3',
        'RNA28S4',
        'RNA28S5',
        'RNA28SN1',
        'RNA28SN2',
        'RNA28SN3',
        'RNA28SN4',
        'RNA28SN5',
        'RNA45S1',
        'RNA45S2',
        'RNA45S3',
        'RNA45S4',
        'RNA45S5',
        'RNA45SN1',
        'RNA45SN2',
        'RNA45SN3',
        'RNA45SN4',
        'RNA45SN5',
        'RNA5-8S1',
        'RNA5-8S2',
        'RNA5-8S3',
        'RNA5-8S4',
        'RNA5-8S5',
        'RNA5-8SN1',
        'RNA5-8SN2',
        'RNA5-8SN3',
        'RNA5-8SN4',
        'RNA5-8SN5',
        'RNA5S1',
        'RNA5S10',
        'RNA5S11',
        'RNA5S12',
        'RNA5S13',
        'RNA5S14',
        'RNA5S15',
        'RNA5S16',
        'RNA5S17',
        'RNA5S2',
        'RNA5S3',
        'RNA5S4',
        'RNA5S5',
        'RNA5S6',
        'RNA5S7',
        'RNA5S8',
        'RNA5S9',
        'RNR1',
        'RNR2',
        'RNR3',
        'RNR4',
        'RNR5'
    ]
    adata.var['gene_group__ribo_rna'] = [
        x in gene_group__ribo_rna for x in adata.var['gene_symbols']
    ]

    # Add in sample metadata.
    # NOTE: it would be more memory efficient to stash this in
    #       unstructured dict-like annotation (adata.uns)


    # Ensure we have experiment_in the final dataframe.
    if 'experiment_id' not in adata.obs.columns:
        adata.obs['experiment_id'] = adata.obs[metadata_key]

    vars_prior_metrics = adata.var_keys()
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=[
            'gene_group__mito_transcript',
            'gene_group__mito_protein',
            'gene_group__ribo_protein',
            'gene_group__ribo_rna'
        ],
        inplace=True
    )


    # Apply cell QC filters.
    adata.obs['cell_passes_hard_filters'] = True
    filters_all_samples = []
    filters_experiment = []
    options.experiment_name = 'all_samples'
    if 'cell_filters' not in params_dict:
        warnings.warn('Found no cell_filters in params_dict.')
    else:
        if options.experiment_name in params_dict['cell_filters'].keys():
            # NOTE: we want this to throw an error if value is not there.
            filters_all_samples = params_dict['cell_filters'][
                options.experiment_name
            ]['value']
        if options.experiment_name in params_dict['cell_filters'].keys():
            filters_experiment = params_dict['cell_filters'][
                options.experiment_name
            ]['value']

    # First record the total number of cells that pass each filter
    # independently i.e., not depenedent on any other filter.
    if len(filters_all_samples) > 0:
        for filter_query in filters_all_samples:
            if filter_query != '':
                try:
                    n_cells_dict['all_cells'][
                        'filter__all_samples {}'.format(filter_query)
                    ] = adata.n_obs - adata.obs.query(filter_query).shape[0]
                except:
                    'specific experiment filter'
    if len(filters_experiment) > 0:
        for filter_query in filters_experiment:
            # NOTE: could add if test here to grab filter_query ==
            # file_cellids_filter or file_cellids_keep
            if filter_query != '':
                n_cells_dict[options.experiment_name][
                    'filter__sample_specific {}'.format(
                        filter_query
                    )
                ] = adata.n_obs - adata.obs.query(filter_query).shape[0]

    # Now apply the filters - first apply the filters for all samples.
    n_cells_start = adata.n_obs
    filter_i = 0
    if len(filters_all_samples) > 0:
        # Run each filter iteratively.
        for filter_query in filters_all_samples:
            if filter_query != '':
                # Drop the cells that are flagged in this query
                cells_to_remove = adata.obs.query(filter_query).index
                adata.obs.loc[cells_to_remove, 'cell_passes_hard_filters'] = False
                # adata = adata[
                #     np.invert(adata.obs.index.isin(cells_to_remove)),
                #     :
                # ]
                if y_n_print:
                    print('[{}] {} "{}": {} dropped {} remain'.format(
                        'all sample cell QC applied',
                        options.experiment_name,
                        filter_query,
                        len(cells_to_remove),
                        adata.obs['cell_passes_hard_filters'].sum()
                    ))
                n_cells_dict[options.experiment_name][
                    'filter__all_samples after_filter_{} {}'.format(
                        filter_i,
                        filter_query
                    )
                ] = adata.obs['cell_passes_hard_filters'].sum()
                filter_i += 1

    # Now apply per sample filters.
    if len(filters_experiment) > 0:
        # Run each filter iteratively.
        for filter_query in filters_experiment:
            if filter_query != '':
                cells_to_remove = adata.obs.query(filter_query).index
                adata.obs.loc[cells_to_remove, 'cell_passes_hard_filters'] = False
                n_cells_dict[options.experiment_name][
                    'filter__sample_specific after_filter_{} {}'.format(
                        filter_i,
                        filter_query
                    )
                ] = adata.obs['cell_passes_hard_filters'].sum()
                filter_i += 1

    # Apply cell downsampling if needed.
    if params_dict['downsample_cells_fraction']['value'] != '':
        n_cells_start = adata.n_obs
        sc.pp.subsample(
            adata,
            fraction=float(
                params_dict['downsample_cells_fraction']['value']
            ),
            copy=False,
            random_state=0
        )
        n_cells_dict[
            options.experiment_name
        ]['downsample_cells_fraction'] = adata.n_obs

    elif params_dict['downsample_cells_n']['value'] != '':
        n_cells_start = adata.n_obs
        sc.pp.subsample(
            adata,
            n_obs=int(params_dict['downsample_cells_n']['value']),
            copy=False,
            random_state=0
        )
        n_cells_dict[
            options.experiment_name
        ]['downsample_cells_n'] = adata.n_obs

    # Apply count downsampling if needed.
    if params_dict['downsample_feature_counts']['value'] != '':
        fraction = params_dict['downsample_feature_counts']['value']
        target_counts_per_cell = adata.obs['total_counts'].apply(
            lambda x: int(x * fraction)
        ).values
        sc.pp.downsample_counts(
            adata,
                counts_per_cell=target_counts_per_cell,
            random_state=0
        )

    # Print the number of cells and genes for this sample.
    n_cells_dict[options.experiment_name]['after_filters'] = adata.obs[
        'cell_passes_hard_filters'
    ].sum()
    
    
    adata.write(
        f'{options.of}.h5ad',
        compression='gzip',
        compression_opts=anndata_compression_opts
    )      
            
if __name__ == '__main__':
    main()