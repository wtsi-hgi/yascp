#!/usr/bin/env python

__date__ = '2020-04-03'
__version__ = '0.0.1'

import argparse
import numpy as np
import pandas as pd
import scanpy as sc


def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Add Sample Metadata to vireo dataset.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-vir', '--vireo',
        action='store',
        dest='vireo',
        required=True,
        help='vireo'
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
        '-esm', '--extra_sample_metadata',
        action='store',
        dest='extra_sample_metadata',
        required=False,
        default='',
        help='extra_sample_metadata'
    )

    options = parser.parse_args()
    vireo = pd.read_csv(options.vireo, sep='\t')
    vireo = vireo.set_index(options.mk)
    extra_sample_metadata = pd.read_csv(options.extra_sample_metadata, sep='\t')
    extra_sample_metadata = extra_sample_metadata.set_index(options.mk)
    for col1 in extra_sample_metadata.columns:
        # print(col1)
        vireo[col1] = extra_sample_metadata[col1]
    vireo.dropna(axis = 1, how = 'all', inplace = True)
    vireo.to_csv('Vireo_metadata.csv',sep='\t')
    print('Done')

if __name__ == '__main__':
    main()