#!/usr/bin/env python

__date__ = '2020-06-26'
__version__ = '0.0.1'
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser(
    description="""
        Basic QC plots.
        """
)

parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s {version}'.format(version=__version__)
)

parser.add_argument(
    '-az', '--az_file',
    action='store',
    dest='az_file',
    required=True,
    help='H5 AnnData file.'
    )

parser.add_argument(
    '-m', '--mapping',
    action='store',
    dest='mapping',
    required=True,
    help='H5 AnnData file.'
    )

parser.add_argument(
    '-of', '--out_file',
    action='store',
    dest='of',
    required=True,
    help='H5 AnnData file.'
    )

options = parser.parse_args()

All_Data = pd.read_csv(options.az_file,compression='gzip',sep='\t')
Mappings = pd.read_csv(options.mapping,sep='\t',index_col=0)

D1 =All_Data
D1['idx1']=D1.index
D1 = All_Data.set_index('predicted.celltype.l2')
for col in Mappings.columns:
    D1[f'{col}_predicted.celltype.l2']=''
    D1[f'{col}_predicted.celltype.l2']=Mappings[col]
D1 =D1.reset_index()
D1 = D1.set_index('idx1')
D1.to_csv(options.of,compression='gzip',sep='\t')    
print('Done')