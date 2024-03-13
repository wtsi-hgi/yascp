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

    # parser.add_argument(
    #     '-pyml', '--params_yaml',
    #     action='store',
    #     dest='pyml',
    #     default='',
    #     help='YAML file containing cell filtering and downsampling (e.g.,\
    #         of total number of cells or reads per cell) parameters.\
    #         If file is not provided, no filtering or downsampling is\
    #         performed.\
    #         (default: %(default)s)'
    # )
    
    parser.add_argument(
        '--downsample_cells_fraction',
        action='store',
        dest='downsample_cells_fraction',
        default='',
        help=''
    )
 
    parser.add_argument(
        '--cell_filters',
        action='store',nargs='+',
        dest='cell_filters',
        default='',
        help=''
    )  
    parser.add_argument(
        '--cell_filters_experiment',
        action='store',nargs='+',
        dest='cell_filters_experiment',
        default='',
        help=''
    ) 

    parser.add_argument(
        '--downsample_feature_counts',
        action='store',
        dest='downsample_feature_counts',
        default='',
        help=''
    )     
 
    parser.add_argument(
        '--downsample_cells_n',
        action='store',
        dest='downsample_cells_n',
        default='',
        help=''
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
    
    parser.add_argument(
        '-drop', '--drop',
        action='store',
        dest='drop',
        default=False,
        help='Key to link metadata to h5addata_file experiment_id column.\
            (default: %(default)s)'
    )    
    
    
    options = parser.parse_args()
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    sc.settings.n_jobs = options.ncpu  # number CPUs

    # Read in the parameters for downsampling and cell filters.
    # if options.pyml == '':
    #     params_dict = {}
    # else:
    #     with open(options.pyml, 'r') as f:
    #         params_dict = yaml.safe_load(f)
    #     params_dict = params_dict['sample_qc']
    params_dict = {}
    params_dict['downsample_cells_n'] = options.downsample_cells_n
    params_dict['downsample_feature_counts'] = options.downsample_feature_counts
    params_dict['downsample_cells_fraction'] = options.downsample_cells_fraction
    params_dict['cell_filters'] = options.cell_filters
    params_dict['cell_filters_experiment'] = options.cell_filters_experiment
    # /software/hgi/pipelines/yascp_versions/yascp_v1.5/conf/qc.conf
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
    if all(params_dict[k] != '' for k in param_filters_check):
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


    # Add in sample metadata.
    # NOTE: it would be more memory efficient to stash this in
    #       unstructured dict-like annotation (adata.uns)


    # Ensure we have experiment_in the final dataframe.
    if 'Exp' not in adata.obs.columns:
        adata.obs['Exp'] = adata.obs[options.mk]


    # Apply cell QC filters.
    adata.obs['cell_passes_hard_filters'] = True
    filters_all_samples = []
    filters_experiment = []
    
    if options.cell_filters =='':
        print('Found no cell_filters in params_dict.')
    else:
        # if options.experiment_name in params_dict['cell_filters']:
            # NOTE: we want this to throw an error if value is not there.
        filters_all_samples = params_dict['cell_filters']
        # if options.experiment_name in options.cell_filters:
        filters_experiment = {}
        for filt1 in  params_dict['cell_filters_experiment'][0].replace('"','').replace("'",'').split(';'):
            sp1 = filt1.split(":")
            filters_experiment[sp1[0]]=sp1[1] + f' and Exp=="{sp1[0]}"'

    # First record the total number of cells that pass each filter
    # independently i.e., not depenedent on any other filter.
    if len(filters_all_samples) > 0:
        for filter_query in filters_all_samples:
            filter_query
            if filter_query != '':
                try:
                    n_cells_dict[options.experiment_name][
                        'filter__all_samples {}'.format(filter_query)
                    ] = adata.n_obs - adata.obs.query(filter_query).shape[0]
                except:
                    'specific experiment filter'
    if len(filters_experiment) > 0:
        for filter_query in filters_experiment:
            # print(filter_query)
            # NOTE: could add if test here to grab filter_query ==
            # file_cellids_filter or file_cellids_keep
            if filters_experiment[filter_query] != '':
                n_cells_dict[filter_query]={}
                n_cells_dict[filter_query][
                    'filter__sample_specific {}'.format(
                        filter_query
                    )
                ] = adata.n_obs - adata.obs.query(filters_experiment[filter_query]).shape[0]

    # Now apply the filters - first apply the filters for all samples.
    n_cells_start = adata.n_obs
    filter_i = 0
    adata.obs.loc[:, 'cell_passes_hard_filters'] = True
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
        
        for filter_query_pre in filters_experiment:
            filter_query = filters_experiment[filter_query_pre]
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
    if options.drop != 'false':
        print('yes, dropping hard qc fails')
        # we either drop of flag the hard filters files
        adata = adata[adata.obs['cell_passes_hard_filters']==True]
    
    # Apply cell downsampling if needed.
    if options.downsample_cells_fraction!= '':
        n_cells_start = adata.n_obs
        sc.pp.subsample(
            adata,
            fraction=float(
                options.downsample_cells_fraction
            ),
            copy=False,
            random_state=0
        )
        n_cells_dict[
            options.experiment_name
        ]['downsample_cells_fraction'] = adata.n_obs

    elif options.downsample_cells_n != '':
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
    if options.downsample_feature_counts != '':
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