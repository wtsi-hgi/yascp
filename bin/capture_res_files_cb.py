#!/usr/bin/env python

__author__ = 'Matiss Ozols'
__date__ = '2022-03-22'
__version__ = '0.0.1'

import os
import glob
import argparse
from os.path import exists
# this code captures the cb files in two levels subdir or subdir/subdir and emits to the subsequent processes.
parser = argparse.ArgumentParser(
    description="""
        Filter and merge 10x data. Save to AnnData object.
        """
)

parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s {version}'.format(version=__version__)
)
# --method ${params.aggregation_method}
parser.add_argument(
    '-res', '--resolution',
    action='store',
    dest='resolution',
    required=True,
    help='',default=False
)

options = parser.parse_args()
resolution = options.resolution

PREFIX1=f"{os.getcwd()}/"
os.mkdir(f"{PREFIX1}captured")
all_vcfs = glob.glob(f'{PREFIX1}tmp1234/cellbender/*/cellbenderFPR_0pt1*/')
all_vcfs12 = glob.glob(f'{PREFIX1}tmp1234/cellbender/*/cellbender_FPR_0pt1*/')
all_vcfs13 = glob.glob(f'{PREFIX1}tmp1234/cellbender/*/cellbender-FPR_0pt1*/')
all_vcfs.extend(all_vcfs12)
all_vcfs.extend(all_vcfs13)
all_vcfs2 = glob.glob(f'{PREFIX1}tmp1234/cellbender/*/*/cellbenderFPR_0pt1*/')
all_vcfs22 = glob.glob(f'{PREFIX1}tmp1234/cellbender/*/*/cellbender_FPR_0pt1*/')
all_vcfs23 = glob.glob(f'{PREFIX1}tmp1234/cellbender/*/*/cellbender-FPR_0pt1*/')
all_vcfs2.extend(all_vcfs22)
all_vcfs2.extend(all_vcfs23)
for vcf1 in all_vcfs2:
    # vcf1=vcf1[:-1]
    vcf1='/'.join(vcf1.split('/')[:-1])
    name = vcf1.split('/')[-3]
    ln2 = f'{PREFIX1}captured/{name}'
    if (exists(ln2)):
        # print('exists')
        _=''
    else:
        # print('doesnt')
        os.mkdir(ln2)
        os.system(f'ln -s {vcf1} {ln2}')
        
for vcf1 in all_vcfs:
    vcf1='/'.join(vcf1.split('/')[:-1])
    name = vcf1.split('/')[-2]
    ln2 = f'{PREFIX1}captured/{name}'
    if (exists(ln2)):
        # print('exists')
        _=''
    else:
        # print('doesnt')
        os.mkdir(ln2)
        os.system(f'ln -s {vcf1} {ln2}')
print('Done')