#!/usr/bin/env python

__date__ = '2023-06-26'
__version__ = '0.0.2'

import pandas as pd
import argparse
import pandas as pd
import re
import glob
# Vireo is not working well if multiple same donor genotypes are provided.
# For this reasone we have to ensure they are removed.
parser = argparse.ArgumentParser(
    description="""
        Joins all the outputs from the GT match in an expected input format for subsequent relatedness estimation.
        """
)

parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s {version}'.format(version=__version__)
)

parser.add_argument(
    '-d', '--donors',
    action='store',
    dest='donors',
    required=True,
    help='donors'
)

parser.add_argument(
    '-b', '--bridge',
    action='store',
    dest='bridge',
    required=True,
    help='bridge'
)

parser.add_argument(
    '-i', '--yascp_input_data',
    action='store',
    dest='yascp_input_data',
    required=True,
    help='yascp_input_data'
)

parser.add_argument(
    '-n', '--name',
    action='store',
    dest='name',
    required=True,
    help='name'
)

options = parser.parse_args()
donors = pd.read_csv(options.donors,sep='\t',header=None)
bridge = pd.read_csv(options.bridge,sep='\t')
input_table=pd.read_csv(options.yascp_input_data,sep='\t')
pool_name = options.name
input_table=pd.read_csv(options.yascp_input_data,sep='\t')
input_expected = input_table.loc[input_table['experiment_id']==pool_name,'donor_vcf_ids'].values[0]
# /lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/fetch/Cardinal_46499_Jan_21_2023/results/yascp_inputs/input.tsv/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/fetch/Cardinal_46499_Jan_21_2023/results/yascp_inputs/input.tsv
newlist = [] # empty list to hold unique elements from the list
duplist = []
mylist = bridge['oragene_id']
for i in mylist:
    if i not in newlist:
        newlist.append(i)
    else:
        duplist.append(i) 

donors['mapping']=''
# donors[0] = donors[0].str.replace("^1_",'')
for don2 in donors[0]:
    # print(don1)
    mappings=[]    
    # don1='1_15002002032964_206862850133_R09C01'
    mappings = bridge[bridge['oragene_id']==don2]['s00046_id']
    
    if len(mappings)>0:
        mapping=mappings.values[0]
        for mapping in mappings:
            print(mapping)
            if mapping in input_expected:
                donors.loc[donors[0]==don2,'mapping']=mapping
                # there are cases where two genotype ids map to same expected donor. In this case we should pick the one expected.
                break
    don1=re.sub('^1_','',don2)        
    mappings = bridge[bridge['oragene_id']==don1.split('_')[0]]['s00046_id']
    if len(mappings)>0:
        mapping=mappings.values[0]
        for mapping in mappings:
            if mapping in input_expected:
                donors.loc[donors[0]==don2,'mapping']=mapping
                # there are cases where two genotype ids map to same expected donor. In this case we should pick the one expected.
                break
        # continue
    try:        
        mappings = bridge[bridge['oragene_id']==don1.split('_')[1]]['s00046_id']
        if len(mappings)>0:
            mapping=mappings.values[0]
            for mapping in mappings:
                if mapping in input_expected:
                    donors.loc[donors[0]==don2,'mapping']=mapping
                    # there are cases where two genotype ids map to same expected donor. In this case we should pick the one expected.
                    break
    except:
        _='error caused by a single name val'
        
donors = donors.sort_values(0)
try:
    g = glob.glob(f'/lustre/scratch123/hgi/projects/cardinal_analysis/freezes/freeze0/qc/*/infered_genotypes/{pool_name}/gt_match_FullGT_Jul2024_final/{pool_name}/GT_replace_stats_{pool_name}_gt_donor_assignments.csv')
    if len(g)>0:
        d2=pd.read_csv(g[0],sep='\t')
        d2=d2[d2['Match Expected']==True]
        all_definite_donors = set(d2['donor_gt original'])
        donors2 = donors[donors[0].isin(all_definite_donors)]
        missing=set(donors['mapping'])-set(donors2['mapping'])
        # missing=set(['30007480246'])
        miss = donors[donors['mapping'].isin(missing)]
        donors = pd.concat([donors2,miss])
except:
    print('this is not available')
    
donors = donors.drop_duplicates(subset=['mapping'])
donors[donors['mapping']==donors['mapping']]
donors = donors[donors['mapping']!='']
donors[0].to_csv('t.tsv',index=False,header=None)
expected = set(input_expected.replace("'",'').split(','))
missing = expected-set(donors['mapping'])
print('Done')


# set(d2['donor_gt original']) - set(donors[0])
# d2=pd.read_csv("/lustre/scratch123/hgi/projects/cardinal_analysis/freezes/freeze0/qc/Cardinal_45327_Jul_18_2022/infered_genotypes/CRD_CMB13016576/gt_match_FullGT_Jul2024_final/CRD_CMB13016576/GT_replace_stats_CRD_CMB13016576_gt_donor_assignments.csv",sep='\t')