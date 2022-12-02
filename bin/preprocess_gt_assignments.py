#!/usr/bin/env python3
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Replace Donor Ids for GT match")
parser.add_argument("-i", "--infile", metavar="donor_id", dest = "in_file",
                    help="", type=str,
                    default=1)

args = parser.parse_args()
Data = pd.read_csv(args.in_file,sep='\t')
Data = Data[Data['Match Expected']==True]
Data['pool_donor'] = Data['donor_query']+'_'+Data['pool']
Data_Mapping = Data[['pool_donor','donor_gt original']]
Data_Mapping.to_csv('genotype_phenotype.tsv',sep='\t',index=False)
print('Done')