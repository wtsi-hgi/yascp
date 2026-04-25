#!/usr/bin/env python

__date__ = '2022-09-20'
__version__ = '0.0.1'

import os
import pandas as pd
import subprocess
import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser(
    description="""
        Takes the required and optional inputs
        """
)
parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s {version}'.format(version=__version__)
)

parser.add_argument(
    '-vcfs', '--vcfs',
    action='store',
    dest='vcfs',
    required=False,
    default=None,
    help='vcfs file to use to convert samples to the mappings'
)

options = parser.parse_args()
List_of_VCfs = options.vcfs


# This code takes all the provided vcf/bcf files and extracts the overlapping positions between them all. 
# It then records it in a bed file, and this can then be used to subset files before merging them all.
# List_of_VCfs = 'sorted_Expected_CRD_CMB13102395.bcf.gz CRD_CMB13102395_headfix_vireo.vcf.gz'
All_VCFs = List_of_VCfs.split(" ")

Overlap_list = set()
for vcf1 in All_VCFs:
    print(vcf1)
    # subprocess.check_output(f"bcftools query -f '%CHROM %POS %REF %ALT\n' {vcf1} |", shell=True, text=True)
    output = subprocess.getoutput(f"bcftools query -f '%CHROM;%POS\n' {vcf1} ")
    All_positions = set(output.split('\n'))

    if len(Overlap_list)==0:
        Overlap_list =All_positions
    else:
        Overlap_list = Overlap_list.intersection(All_positions)
Bed_File = pd.DataFrame(Overlap_list,columns=['chr_pos'])
Bed_File['chr']=Bed_File['chr_pos'].str.split(';').str[0]
Bed_File['pos']=Bed_File['chr_pos'].str.split(';').str[1].astype(int)-1
Bed_File['pos2']=Bed_File['chr_pos'].str.split(';').str[1].astype(int)
Bed_File_record = Bed_File[['chr','pos','pos2']]
# Bed_File_record = Bed_File['chr']+':'+Bed_File['pos'].astype(str)
Bed_File_record.to_csv('Bed_File_record.bed',header=False,index=False,sep='\t')
# bcftools view -R Bed_File_record.bed -Ov -o  test.vcf sorted_Expected_CRD_CMB13102395.bcf.gz
print('Done')

