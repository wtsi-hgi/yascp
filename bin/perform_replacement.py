#!/usr/bin/env python

__date__ = '2020-06-30'
__version__ = '0.0.1'


import pandas as pd
import argparse
import numpy as np

from scipy import stats

parser = argparse.ArgumentParser(description="Replace Donor Ids for GT match")
parser.add_argument("-id", "--donor_id", metavar="donor_id", dest = "donor_id",
                    help="", type=str,
                    default=1)
parser.add_argument("-gp", "--genotype_phenotype_mapping", metavar="genotype_phenotype_mapping", dest = "genotype_phenotype_mapping",
                    help="", type=str,
                    default=None)

parser.add_argument("-a", "--assignemts", metavar="assignemts", dest = "assignemts",
                    help="", type=str,
                    default=None)

parser.add_argument("-gt", "--gt", metavar="gt", dest = "gt",
                    help="", type=str,
                    default=None)


args = parser.parse_args()
if (args.genotype_phenotype_mapping):
    gp_ma = pd.read_csv(args.genotype_phenotype_mapping,sep='\t',index_col=0)
    gp_ma.index = gp_ma.index.astype(str)
    # gt = pd.read_csv(args.gt)
    # gt.iloc[:,0].str.split('\t')

    # here we will perform the substitution and overwrite the GT file
    assignemts = pd.read_csv(args.assignemts,sep=',')
    assignemts=assignemts.set_index('donor_gt')
    assignemts['donor_gt']=''
    for ix in assignemts.index:
        print(ix)
        try:
            replacement = gp_ma.loc[ix,'donor_gt'].values[0]
        except:
            # try:
            split2 = ix.split('_')[0]
            replacement = gp_ma.loc[split2].values[0]
                
            #     print(replacement)
            # except:
            #     replacement = f'No mapping: {ix}'
        # print(replacement)
        assignemts.loc[ix,'donor_gt']=replacement

        # sp1 = gt.iloc[:,0].str.split('\t').str[2].replace(ix, replacement)
        # gt.iloc[:,0] = gt.iloc[:,0].str.split('\t').str[:2].str.join('\t') +'\t'+sp1+'\t'+gt.iloc[:,0].str.split('\t').str[3:].str.join('\t') 
    assignemts = assignemts.set_index('donor_gt')
    assignemts = assignemts.reset_index()
    assignemts = assignemts[['donor_query','donor_gt','score0','score1','score_n','n','mean','sd','z0','z1']]
    assignemts.to_csv(f'pheno_{args.assignemts}',sep=',',index=False)

    print('Done')
else:
    gp_ma = pd.DataFrame()

