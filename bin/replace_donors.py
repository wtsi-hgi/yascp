#!/usr/bin/env python3
import pandas as pd
import argparse
import numpy as np

from scipy import stats

parser = argparse.ArgumentParser(description="Replace Donor Ids for GT match")
parser.add_argument("-id", "--donor_id", metavar="donor_id", dest = "donor_id",
                    help="", type=str,
                    default=1)
args = parser.parse_args()
input_id =args.donor_id
# input_id = 'ELGH_VAL11650907'
GT_Assignments = pd.read_csv(f'{input_id}_assignments.csv',sep=',',index_col=0)

GT_Assignments__sample_summary_txt = pd.read_csv(f'{input_id}.sample_summary.txt',sep='\t',header=None)
GT_Assignments__exp_sample_summary_txt = pd.read_csv(f'{input_id}__exp.sample_summary.txt',sep='\t',header=None)
donor_ids=pd.read_csv(f'donor_ids.tsv',sep='\t')
vcf = pd.read_csv(f'GT_donors.vireo.vcf',sep='XXXXXXXdummy')

slipt1 = vcf[:1].values[0,0].split('\t')
# detect the outlier z values and make this as unassigned.

from scipy.stats import iqr
iqr = iqr(GT_Assignments['z_score'])
q3 = np.percentile(GT_Assignments['z_score'], 75) 
q1 = np.percentile(GT_Assignments['z_score'], 25) 
outlier_thresh = q3+1.5*iqr
outlier_thresh_neg = q1-1.5*iqr

Unassigned = list(GT_Assignments[GT_Assignments['z_score']<outlier_thresh_neg].index)


for ix in GT_Assignments.index:
    print(ix)
    #Here should add a a filter to estimate whether it is a good match.
    if ix in Unassigned:
       replacement = 'Unassigned'
    else: 
        replacement = GT_Assignments.loc[ix,'donor_gt']
    GT_Assignments__exp_sample_summary_txt.loc[:,0] = GT_Assignments__exp_sample_summary_txt.loc[:,0].replace(f'{input_id}__{ix}', f'{input_id}__{replacement}')
    GT_Assignments__sample_summary_txt.loc[:,1]=GT_Assignments__sample_summary_txt.loc[:,1].replace(ix, replacement)
    donor_ids['donor_id']=donor_ids['donor_id'].replace(ix, replacement)
    donor_ids['best_singlet']=donor_ids['best_singlet'].replace(ix, replacement)
    sp1 = donor_ids['best_doublet'].str.split(',').str[0].replace(ix, replacement)
    sp2 = donor_ids['best_doublet'].str.split(',').str[1].replace(ix, replacement)
    donor_ids['best_doublet'] =sp1+','+sp2
    slipt1=list(map(lambda x: x.replace(ix, replacement), slipt1))

slipt1='\t'.join(slipt1)
vcf[:1].values[0,0]=slipt1
# Now that we have loaded them replace all the files with the correct donor ids.

GT_Assignments__sample_summary_txt.to_csv(f'GT_replace_{input_id}.sample_summary.txt',sep='\t',index=False)
GT_Assignments__exp_sample_summary_txt.to_csv(f'GT_replace_{input_id}__exp.sample_summary.txt',sep='\t',index=False)
donor_ids.to_csv(f'GT_replace_donor_ids.tsv',sep='\t',index=False)
vcf.to_csv(f'GT_replace_GT_donors.vireo.vcf',sep=',',index=False)

print("Done")