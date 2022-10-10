#!/usr/bin/env python3
import pandas as pd
import argparse


parser = argparse.ArgumentParser(description="Replace Donor Ids for GT match")
parser.add_argument("-i", "--in_file", metavar="donor_id", dest = "in_file",
                    help="", type=str,
                    default=1)
parser.add_argument("-gp", "--genotype_phenotype_mapping", metavar="genotype_phenotype_mapping", dest = "genotype_phenotype_mapping",
                    help="", type=str,
                    default=None)

parser.add_argument("-if", "--input_file", metavar="input_file", dest = "input_file",
                    help="", type=str,
                    default=None)
                    
args = parser.parse_args()
if (args.genotype_phenotype_mapping):
    gp_ma = pd.read_csv(args.genotype_phenotype_mapping,sep='\t',index_col=0)
    gp_ma.index = gp_ma.index.astype(str)
else:
    gp_ma = pd.DataFrame()
input_id =args.in_file

try:
    GT_Assignments = pd.read_csv(f'{input_id}',sep=',')
except:
    GT_Assignments = pd.read_csv(f'{input_id}',sep='\t')

for i,row1 in GT_Assignments.iterrows():
    # print(row1)
    math=row1['donor_gt original']
    # math = 'celline_THP1'
    expected = str(row1['Emergency_ids expected'])+','+str(row1['Good_ids expected'])
    expected= expected.split(',')
    DF2 = pd.DataFrame(expected, columns=['col1'])
    try:
        repl = str(gp_ma.loc[math.split('_')[0]][0])
        GT_Assignments.loc[i,'donor_gt']=repl
        if (len(DF2[DF2.col1.str.contains(repl)])>0):
            GT_Assignments.loc[i,'Match Expected']=True
    except:
        try:
            repl = str(gp_ma.loc[math][0])
            GT_Assignments.loc[i,'donor_gt']=repl
            if (len(DF2[DF2.col1.str.contains(repl)])>0):
                GT_Assignments.loc[i,'Match Expected']=True
        except:
            print('no need to replace')
GT_Assignments.to_csv(input_id,sep='\t',index=False)