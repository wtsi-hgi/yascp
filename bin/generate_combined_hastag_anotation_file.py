#!/usr/bin/env python

__date__ = '2021-11-15'
__version__ = '0.0.1'

import argparse
from distutils.version import LooseVersion
import os
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
import pandas as pd
import scanpy


def combine_reports(all_alternitive,mode):
    all_indexes_full=set({})
    for d1 in all_alternitive:
        pool = d1.split('__')[0]
        if d1 in ('fake_file.fq', 'fake_file1.fq', 'fake_file2.fq'):
            Dataset = pd.DataFrame()
        else:
            Dataset = pd.read_csv(d1,sep='\t',index_col=0)
            
            if(len(Dataset.columns)==0):
                Dataset = pd.read_csv(d1,sep=',',index_col=0)
            Dataset.index = Dataset.index +'-'+pool
        Dataset=Dataset.add_prefix(mode)
        
        all_indexes = set(Dataset.index)
        all_indexes_full = all_indexes_full.union(all_indexes)
    Data_All_alt=pd.DataFrame(index=list(set(all_indexes_full)))    
    for d1 in all_alternitive:
        pool = d1.split('__')[0]
        if d1 in ('fake_file.fq', 'fake_file1.fq', 'fake_file2.fq'):
            Dataset = pd.DataFrame()
        else:
            Dataset = pd.read_csv(d1,sep='\t',index_col=0)
            if(len(Dataset.columns)==0):
                Dataset = pd.read_csv(d1,sep=',',index_col=0)
            Dataset.index = Dataset.index +'-'+pool
        Dataset=Dataset.add_prefix(mode)
        for col1 in Dataset.columns:
            try:
                _ = Data_All_alt[col1]
            except:
                Data_All_alt[col1]=''
            Data_All_alt.loc[Dataset.index,col1] = Dataset[col1]
    return Data_All_alt

def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Merges all the hastag information in one file for hastag and Celltypist.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-az', '--all_hastag_files',
        action='store',
        dest='all_hastag_files',
        required=True,
        help='List of csv-delimited files of hastag assigned by hastag.'
    )


    options = parser.parse_args()
    
    Data_All=pd.DataFrame()
    # Read hastag files from a TSV file using pandas
    hastag_df = pd.read_csv(options.all_hastag_files, header=None, names=['file_path'])
    hastag_files = hastag_df['file_path'].tolist()
    Data_All = combine_reports(hastag_files,'hastag:')
    Donor_Exp = Data_All.index.map(lambda x: '-'.join(x.split('-')[2:]))
    Donor = Donor_Exp.str.split('__').str[-1]
    Exp = Donor_Exp.str.split('__').str[0]
    Data_All['Donor'] =Donor
    Data_All['Exp'] =Exp
    Data_All.to_csv('All_hastag_Assignments.tsv',sep='\t')


if __name__ == '__main__':
    main()
