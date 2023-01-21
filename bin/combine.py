#!/usr/bin/env python3


__date__ = '2021-04-11'
__version__ = '0.0.1'

import argparse
import pandas as pd



def main(files,outfile,seed):
    dataframe=pd.DataFrame()
    for file in files.split(" "):
        if seed:
            cols = seed.split(',')
            Data= pd.read_csv(file,sep='\t',header=None, names = cols)
        else:
            Data= pd.read_csv(file,sep='\t')
        dataframe = pd.concat([dataframe,Data], sort=False)
    cols = dataframe.columns
    print(cols)
    dataframe.to_csv(outfile,sep='\t',index=False,columns=dataframe.columns)

if __name__ == '__main__':


    parser = argparse.ArgumentParser(
        description="""
            Combines all the Raw cellranger metadata to be passed in the h5ad
            """
    )

    parser.add_argument(
        '-d',
        action='store',
        dest='paths',
        required=True,
        help='Metadata Paths')

    parser.add_argument(
        '-seed',
        action='store',
        dest='seed',
        required=False,
        default=None,
        help='Metadata Paths')
    
    parser.add_argument(
        '-o',
        action='store',
        dest='outfile',
        required=True,
        help='Metadata Paths')
        
    options = parser.parse_args()
    outfile = options.outfile
    files = options.paths
    seed = options.seed
    main(files,outfile,seed)
