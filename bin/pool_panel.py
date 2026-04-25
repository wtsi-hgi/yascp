#!/usr/bin/env python
import pandas as pd
import os
import argparse
import string
import random

# Purpose of this code is to emit the genotypes so that we do not replicate the data. All the genotypes with unique compositions will be processed.
# If more than one pool has the same composition we create a random hex and mapping to tos genotype for each pool is indicated in a tsv file.

parser = argparse.ArgumentParser(description="Replace Donor Ids for GT match")
parser.add_argument("-ids", "--pool_ids", metavar="pool_ids", dest = "pool_ids",
                    help="", type=str,
                    default='')
parser.add_argument("-vcf", "--vcf", metavar="vcf", dest = "vcf",
                    help="", type=str,
                    default='')
parser.add_argument("-mode", "--mode", metavar="mode", dest = "mode",
                    help="", type=str,
                    default='')

args = parser.parse_args()
pool_ids = args.pool_ids
vcfs = args.vcf

def id_generator(size=6, chars=string.ascii_uppercase + string.digits):
    return ''.join(random.choice(chars) for _ in range(size))

# pool_ids = '_pool1::pool2::pool3::pool4::pool5::pool6::pool7::pool8::pool9::pool10::pool11::pool12::pool13::pool14::pool15::pool16::pool17::pool18::pool19::pool20::pool21::pool22::pool23::pool24::pool25::pool26::pool27::pool28::pool29::pool30::pool31::pool32::pool33::pool34::pool35::pool36::pool37::pool38::pool39::pool41::pool42::pool43::pool44::pool45::pool46::pool47::pool48::pool49::pool50::pool51::pool52::pool53::pool54::pool55::pool56::pool57::pool58::pool59::pool60::pool61::pool62::pool63::pool64::pool65::pool67::pool68::pool69::pool70::pool71::pool72::pool73::pool74::pool75::pool76::pool77'
# vcfs = 'Study_Merge_AllExpectedGT_4ZJ8410RF_out.vcf.gz'
if (pool_ids[0]=='_'):
    pool_ids = pool_ids[1:]
pool_ids=pool_ids.split('::')
# We define a random index pool name if more than onoe pool contains the required genotype, otherwise we use rand index
if len(pool_ids)>1:
    Genotype_folder = f"Genotype___{args.mode}_{id_generator(9)}"
else:
    Genotype_folder = f"Genotype___{args.mode}_{pool_ids[0]}"

# Make the folder and Symlink the genotype in tnere 
# os.chdir('/lustre/scratch125/humgen/teams/hgi/mo11/oneK1k/work/a0/dbf245e3d2adddbdf817dc7b7b50c1')
cwd1 = os.getcwd()
os.mkdir(Genotype_folder)
os.system(f'cd ./{Genotype_folder} && ln -s ../{vcfs} && ln -s ../{vcfs}.csi' )
    
Dataset_Pipeline=[]
Dataset_User=[]
for pid in pool_ids:
    Dataset_Pipeline.append({'Pool_id':pid,'vcf':f'{cwd1}/{Genotype_folder}/{vcfs}','vcf_csi':f'{cwd1}/{Genotype_folder}/{vcfs}.csi'})
    Dataset_User.append({'Pool_id':pid,'vcf':f'{Genotype_folder}/{vcfs}'})
Data_Pipeline = pd.DataFrame(Dataset_Pipeline)
Data_User = pd.DataFrame(Dataset_User)

Data_Pipeline.to_csv(f'Data_Pipeline___{Genotype_folder}.tsv',sep='\t',index=False)
Data_User.to_csv(f'Data_User___{Genotype_folder}.tsv',sep='\t',index=False)

print('Done')

