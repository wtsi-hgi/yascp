#!/usr/bin/env python
import os
import glob
import argparse
from os.path import exists
# this code captures the cb files in two levels subdir or subdir/subdir and emits to the subsequent processes.
# parser = argparse.ArgumentParser(
#     description="""
#         Filter and merge 10x data. Save to AnnData object.
#         """
# )

PREFIX1=f"{os.getcwd()}/"
os.mkdir(f"{PREFIX1}captured")
all_vcfs = glob.glob(f'{PREFIX1}tmp1234/cellbender/*/cellbender-FPR_0pt1*')
all_vcfs2 = glob.glob(f'{PREFIX1}tmp1234/cellbender/*/*/cellbender-FPR_0pt1*')
for vcf1 in all_vcfs2:
    # print(vcf1)
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
    # print(vcf1)
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