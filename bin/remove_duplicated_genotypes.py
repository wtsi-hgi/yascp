#!/usr/bin/env python

__date__ = '2020-03-13'
__version__ = '0.0.2'

import pandas as pd
import argparse
import pandas as pd
import re
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
options = parser.parse_args()

donors = pd.read_csv(options.donors,sep='\t',header=None)
bridge = pd.read_csv(options.bridge,sep='\t')

donors['mapping']=''
for don1 in donors[0]:
    print(don1)
    mappings = bridge[bridge['oragene_id']==don1]['s00046_id']
    if len(mappings)>0:
        mapping=mappings.values[0]
        donors.loc[donors[0]==don1,'mapping']=mapping
        continue
    mappings = bridge[bridge['oragene_id']==don1.split('_')[0]]['s00046_id']
    if len(mappings)>0:
        mapping=mappings.values[0]
        donors.loc[donors[0]==don1,'mapping']=mapping
        continue
        
donors = donors.sort_values(0)
donors = donors.drop_duplicates(subset=['mapping'])
donors[0].to_csv('t.tsv',index=False,header=None)
print('Done')

