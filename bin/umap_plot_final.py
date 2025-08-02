#!/usr/bin/env python

__date__ = '2021-22-10'
__version__ = '0.0.2'

import argparse
import os
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib import gridspec

# Nice large palette.
COLORS_LARGE_PALLETE = [
    '#0F4A9C', '#3F84AA', '#C9EBFB', '#8DB5CE', '#C594BF', '#DFCDE4',
    '#B51D8D', '#6f347a', '#683612', '#B3793B', '#357A6F', '#989898',
    '#CE778D', '#7F6874', '#E09D37', '#FACB12', '#2B6823', '#A0CC47',
    '#77783C', '#EF4E22', '#AF1F26'
]


# This function is based off of scanpy:
# https://github.com/theislab/scanpy/blob/master/scanpy/plotting/_tools/scatterplots.py
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


def save_plot(
    adata,
    dict__umap_dim_and_params,
    out_file_base,
    color_var,
    colors_quantitative=True,
    colors_large_palette=COLORS_LARGE_PALLETE,
    drop_legend=25
):
    """Save a plot."""
    n_params = len(dict__umap_dim_and_params[
        list(dict__umap_dim_and_params.keys())[0]
    ])
    # fig_size = 0.125*n_params
    # if fig_size>2^16-10:
    #     fig_size=2^15
    # fig, grid = panel_grid(
    #     hspace=fig_size,
    #     wspace=None,
    #     ncols=4,
    #     num_panels=len(dict__umap_dim_and_params)
    # )
    # i__ax = 0
    fig, grid = panel_grid(
        hspace=0.125*n_params,
        wspace=None,
        ncols=4,
        num_panels=len(dict__umap_dim_and_params)
    )
    i__ax = 0
    # try:
    #     import pandas as pd
    #     import json
    #     PCAs =pd.DataFrame(adata.obsm['X_pca']).set_index(adata.obs.index).iloc[:,1:20]
    #     UMAPs =pd.DataFrame(adata.obsm['X_umap']).set_index(adata.obs.index)
    #     params = pd.DataFrame(adata.obs)
    #     data_dict={}
    #     data_dict["PCAs"]=PCAs.to_dict()
    #     data_dict["UMAPs"]=UMAPs.to_dict()
    #     data_dict["params"]=params.to_dict()
    #     with open(f"Data.json", 'w') as json_file:
    #         json.dump(data_dict, json_file)
    # except:
    #     print('no metadata.s')
    fig.suptitle(color_var, fontsize=16)
    for i__umap in dict__umap_dim_and_params:

        # Set the facet title.
        plt_title = ''
        for parameter, value in dict__umap_dim_and_params[i__umap].items():
            plt_title = '{}{}={}\n'.format(
                plt_title,
                parameter,
                value
            )
        plt_title = plt_title.rstrip()

        # Get the proper axis for this plot.
        ax = plt.subplot(grid[i__ax])

        # We could avoid this line by setting basis=i__umap, but then not
        # consistent axis labels (e.g., X_umap__n_neighbors_151).
        adata.obsm['X_umap'] = adata.obsm[i__umap]

        # If not colors_quantitative, then boot up categorical plots.
        legend_loc = 'right margin'
        color_palette = 'tab20'
        if colors_quantitative is False:
            # Cast to category - required for booleans.
            adata.obs[color_var] = adata.obs[color_var].astype('category')
            n_categories = len(adata.obs[color_var].cat.categories)
            color_palette = None
            if n_categories <= len(plt.get_cmap('Dark2').colors):
                color_palette = 'Dark2'
            elif n_categories <= len(colors_large_palette):
                color_palette = colors_large_palette
            if drop_legend >= 0 and n_categories >= drop_legend:
                legend_loc = None
        else:
            # Make sure that the numeric value is actually numeric
            adata.obs[color_var] = adata.obs[color_var].astype('double')

        # For some reason, there ics an error if return_fig = False, but
        # not when show = False.
        # sc.pl.embedding(
        #     basis='X_umap',
        # Azimuth:predicted.celltype.l2.score
        
        if color_var != 'embedding_density':
            sc.pl.umap(
                adata=adata,
                color=color_var,
                palette=color_palette,
                alpha=0.4,
                title=plt_title,
                legend_loc=legend_loc,
                ax=ax,
                show=False
            )
        else:
            # NOTE: i__umap looks something like X_umap__n_neighbors=15...
            adata.obs['umap_density'] = adata.obs[
                '{}__density'.format(i__umap).replace('X_', '')
            ]
            adata.uns['umap_density_params'] = adata.uns[
                '{}__density_params'.format(i__umap).replace('X_', '')
            ]
            sc.pl.embedding_density(
                adata=adata,
                basis='umap',
                alpha=0.4,
                title=plt_title,
                ax=ax,
                show=False
            )
            del adata.obs['umap_density']
            del adata.uns['umap_density_params']

        del adata.obsm['X_umap']
        i__ax += 1

    fig.savefig(
        'umap-{}.png'.format(color_var),
        #dpi=300,
        bbox_inches='tight'
    )


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Read AnnData object and PCs file. Generates UMAP.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-h5', '--h5_anndata',
        action='store',
        dest='h5',
        required=True,
        help='H5 AnnData file.'
    )

    # parser.add_argument(
    #     '-h5_QC', '--h5_anndata_QC',
    #     action='store',
    #     dest='h5_QC',
    #     required=True,
    #     help='H5 AnnData file.'
    # )

    parser.add_argument(
        '-cq', '--colors_quantitative',
        action='store',
        dest='cq',
        default='',
        help='Comma seperated list of quantitative variable names for colors.\
            (default: "")'
    )

    parser.add_argument(
        '-cc', '--colors_categorical',
        action='store',
        dest='cc',
        default='',
        help='Comma seperated list of categorical variable names for colors.\
            (default: "")'
    )

    parser.add_argument(
        '-dln', '--drop_legend_n',
        action='store',
        dest='drop_legend',
        default=-1,
        type=int,
        help='Drop the legend for categorical colors with >= drop_legend_n\
            categories. If drop_legend_n < 0, then no legend drops are\
            performed.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-ncpu', '--number_cpu',
        action='store',
        dest='ncpu',
        default=4,
        type=int,
        help='Number of CPUs to use.\
            (default: %(default)s)'
    )

    parser.add_argument(
        '-of', '--output_file',
        action='store',
        dest='of',
        default='',
        help='Basename of output files, assuming output in current working \
            directory.\
            (default: <h5_anndata>-<tsv_pcs>-umap)'
    )

    options = parser.parse_args()

    # Scanpy settings
    sc.settings.figdir = os.getcwd()  # figure output directory to match base.
    sc.settings.n_jobs = options.ncpu  # number CPUs
    # sc.settings.max_memory = 500  # in Gb
    #sc.set_figure_params(dpi_save=300)

    # Load the AnnData file.
    adata = sc.read_h5ad(filename=options.h5)
    # adata_QC = sc.read_h5ad(filename=options.h5_QC)
    # adata_QC = sc.read_h5ad(filename='/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/Franke_with_genotypes/work/c5/2ead02a35b945b2fded83b6402ffee/adata_celltypes.h5ad')
    try:
        del adata.layers['counts']
    except:
        _='doesnt exist'
 
    try:
        del adata.layers['log1p_cp10k']
    except:
        _='doesnt exist'
    
    # adata.obs['predicted.celltype.l2']=adata_QC.obs['predicted.celltype.l2']
    # adata.obs['predicted.celltype.l2.score']=adata_QC.obs['predicted.celltype.l2.score']
    # adata.obs['mapping.score']=adata_QC.obs['mapping.score']
    # adata.obs['prob_doublet']=adata_QC.obs['prob_doublet']
    
    # Get the out file base.
    out_file_base = options.of
    if out_file_base == '':
        out_file_base = '{}-{}-umap'.format(
            os.path.basename(options.h5.rstrip('h5ad').rstrip('.'))
        )

    import re

    colors_quantitative = []
    if options.cq != '':
        for pattern in options.cq.split(','):
            regex = re.compile(pattern)
            matches = [col for col in adata.obs.columns if regex.search(col)]
            colors_quantitative.extend(matches)

    colors_categorical = []
    if options.cc != '':
        for pattern in options.cc.split(','):
            regex = re.compile(pattern)
            matches = [col for col in adata.obs.columns if regex.search(col)]
            colors_categorical.extend(matches)

    # Optional: remove duplicates
    colors_quantitative = list(set(colors_quantitative))
    colors_categorical = list(set(colors_categorical))
    
    if len(colors_quantitative) == 0 and len(colors_categorical) == 0:
        raise Exception('Specify a color value.')

    # Add colors_large_palette to adata.uns.
    # adata.uns["annotation_colors"] = COLORS_LARGE_PALLETE

    # Get a list of umap keys
    list__umap_keys = [key for key in adata.obsm if 'X_umap__' in key]
    # Make a dictionary where key = obsm key and value = dictionary of
    # parameters used to generate that obsm.
    dict__umap_dim_and_params = {}
    for key in list__umap_keys:
        dict__umap_dim_and_params[key] = {}
        # Loop over parameters and filter out those we do not care about.
        for param in adata.uns[key.replace('X_umap', 'umap_params')]:
            if param not in set(['a', 'b', 'density']):
                dict__umap_dim_and_params[key][param] = adata.uns[
                    key.replace('X_umap', 'umap_params')
                ][param]

    # NOTE: If the color var is a gene, you should color by ln(CPM+1).
    #       By default these sc.pl.umap uses the .raw attribute of AnnData
    #       if present which is assumed to be ln(CPM+1).
    # colors_quantitative.append('prob_doublet')
    # colors_quantitative.append('pct_counts_gene_group__ribo_rna')
    # colors_quantitative.append('predicted.celltype.l2.score')
    # colors_quantitative.append('mapping.score')
    # colors_categorical.append('predicted.celltype.l2')

    # For each color to plot, loop over the different iterations.

    
    for color_var in colors_quantitative:
        try:
            save_plot(
                adata=adata,
                dict__umap_dim_and_params=dict__umap_dim_and_params,
                out_file_base=out_file_base,
                color_var=color_var,
                colors_quantitative=True,
                drop_legend=options.drop_legend
            )
        except:
            print(f'{color_var} -- cant be plotted')
    for color_var in colors_categorical:
        try:
            save_plot(
                adata=adata,
                dict__umap_dim_and_params=dict__umap_dim_and_params,
                out_file_base=out_file_base,
                color_var=color_var,
                colors_quantitative=False,
                drop_legend=options.drop_legend
            )
        except:
            print(f'{color_var} -- cant be plotted')

    # adata.write(
    #     '{}.h5ad'.format('test'),
    #     compression='gzip'
    # )
    # In some ocassions, you might still observe disconnected clusters and
    # similar connectivity violations. They can usually be remedied by running:
    # sc.tl.paga(adata)
    # From below, remove `plot=False` if you want to see the coarse-grained
    # graph
    # sc.pl.paga(adata, plot=False)
    # sc.tl.umap(adata, init_pos='paga')


if __name__ == '__main__':
    main()
