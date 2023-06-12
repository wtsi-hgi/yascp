#!/usr/bin/env python3
import pandas as pd
import argparse
__date__ = '2021-11-04'
__version__ = '0.0.1'

parser = argparse.ArgumentParser(
    description="""
        Subsetting to 80% variants
        """
)

parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s {version}'.format(version=__version__)
)

parser.add_argument(
    '-random_state', '--random_state',
    action='store',
    dest='random_state',
    default=1,type=int
)

parser.add_argument(
    '-vcf', '--vcf',
    action='store',
    dest='vcf',
    default=None,
)

options = parser.parse_args()

VCF = pd.read_csv(options.vcf,sep='\t',comment='#',header=None)
VCF2 = VCF.sample(frac =.80,random_state=options.random_state)
VCF2=VCF2.sort_index()
VCF2.to_csv('random_variants.tsv',sep='\t',index=False,header=False)
print('Done')