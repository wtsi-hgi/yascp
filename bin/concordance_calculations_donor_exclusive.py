#!/usr/bin/env python3

__date__ = '2023-04-14'
__version__ = '0.0.1'
# python -m debugpy --listen 0.0.0.0:5678 --wait-for-client /lustre/scratch125/humgen/teams/hgi/mo11/oneK1k/extra0/yascp/bin/concordance_calculations_donor_exclusive_work.py --cpus 6 --cell_vcf cellSNP.cells.vcf.gz --donor_assignments stats_pool12_gt_donor_assignments.csv --gt_match_vcf Study_Merge_GTMatchedSubset_EUY8DDDZD_out.vcf.gz --expected_vcf Study_Merge_AllExpectedGT_SYIDTL7VN_out.vcf.gz --cell_assignments GT_replace_donor_ids_true.tsv
import argparse
import sys
import importlib.util
import random
import pickle 
import pandas as pd
import gzip
import numpy as np
import time
import multiprocessing as mp
from multiprocessing import Lock
import logging
import os

class VCF_Loader:
    
    def __init__(self, vcf_file, biallelic_only=True,
                        sparse=False, format_list=['GT']):
        self.vcf_file = vcf_file
        self.load_sample = True
        self.biallelic_only = biallelic_only
        self.sparse = sparse
        self.record_dict={}
        self.reset()
        self.format_list = format_list
        self.exclusive_donor_variants = {}
        self.curently_pushing =[] #this is a lock value to check if rhe curent field is updated so to avaid the race for update
        self.last_count=-1
        self.reset_c()
    
    def reset_c(self):
        self.record_times=0
        
    def reset(self):
        self.exclusive_donor_variants ={}
                
    def myfunc(self):
        print(f"Hello my name is {self.biallelic_only}" )
        
    def load_sample_mp(self,line,obs_ids,count,format_list):
        list_val = line.rstrip().split("\t") #[:5] #:8
        idx = find(list_val[8].split(':'),format_list[0])[0]
        if len(list_val[3]) > 1 or len(list_val[4]) > 1:
            # CURRENTLY DEALS ONLY WITH BIALELIC
            print(f'{idx} var not bialelic')
        else:
            list_val2 = list_val[9:]
            obs = pd.DataFrame(obs_ids)
            lv = pd.DataFrame(list_val2)
            lv_proc =lv[0].str.split(':').str[idx]
            gt_exists = lv_proc[lv_proc != '.']
            idx2 = gt_exists.index
            obs_with_gt = obs.loc[idx2.values]
            obs_with_gt = list(obs_with_gt[0].values)
            list_val_with_gt = lv.loc[idx2.values]
            list_val_with_gt = list(list_val_with_gt[0].values)
            random.seed(count)
            c = list(zip(obs_with_gt, list_val_with_gt))
            random.shuffle(c)
            obs_with_gt, list_val_with_gt = zip(*c)
            # self.append_results([obs_with_gt,list_val_with_gt,idx,list_val,count])

        return [obs_with_gt,list_val_with_gt,idx,list_val,count]

    def set_results(self,to_set,id):
        # Recod to disk to save the loading mmeory time.
        with open(f'tmp_{id}.pkl', 'wb') as f:
            pickle.dump(to_set, f)
        self.record_dict[id]=f'tmp_{id}.pkl'
    
    def append_results(self,result):
        # exclusive_donor_variants
        obs_with_gt= result[0]
        list_val_with_gt= result[1]
        idx = result[2]
        list_val = result[3]
        count = result[4]
        # print(count)
        # print(self.record_times)
        # if
        # print(self.last_count)
        count11=0
        # r = random.random()
        # Issue is that this slows down after number of entries is recorded. So recoding takes longer and longer.
        # every 500 itterations we push the data to a dictionary, later we combine these together.
        if (count % 300 == 0):
            print(f'recording and resetting memory {count}')
            # self.record_dict[count]=self.exclusive_donor_variants
            self.set_results(self.exclusive_donor_variants,count)
            self.reset()  
            self.reset_c()        
        
        for ob_id in obs_with_gt:
            donor_loc_in_list = count11
            alleles = list_val_with_gt[donor_loc_in_list].split(':')[idx]
            if alleles!='.':
                ids = "_".join([list_val[x] for x in [0, 1, 3, 4]])
                donor_var = f"{ids}_{alleles}"
                while ob_id in self.curently_pushing:
                    time.sleep(r*0.01)
                self.curently_pushing.append(ob_id)           
                try:
                    self.exclusive_donor_variants[ob_id].add(donor_var)
                    self.record_times=self.record_times+1
                except:
                    self.exclusive_donor_variants[ob_id]=set()
                    self.exclusive_donor_variants[ob_id].add(donor_var)
                    self.record_times=self.record_times+1
                self.curently_pushing.remove(ob_id)
                # self.exclusive_donor_variants['CTGAAACGTAAGTTCC-1']
            count11+=1 

    def combine_written_files(self):
        to_export = self.exclusive_donor_variants
        for val1 in self.record_dict.values():
            # here remove the int files.
            print(f"merging temp file: {val1}")
            with open(val1, 'rb') as f:
                loaded_dict = pickle.load(f)
                for k1 in loaded_dict.keys():
                    try:
                        to_export[k1]=to_export[k1].union(loaded_dict[k1])
                    except:
                        to_export[k1]=set()
                        to_export[k1]=to_export[k1].union(loaded_dict[k1])
            os.remove(val1)
        return to_export
    
    
    def load_VCF_batch_paralel(self):
        """
        Load whole VCF file by utilising multiple cores to speed up loading of large cell files
        -------------------
        Initially designed to load VCF from cellSNP output, requiring 
        1) all variants have the same format list;
        2) a line starting with "#CHROM", with sample ids.
        If these two requirements are satisfied, this function also supports general
        VCF files, e.g., genotype for multiple samples.

        Note, it may take a large memory, please filter the VCF with bcftools first.
        """
        
        vcf_file = self.vcf_file
        biallelic_only = self.biallelic_only
        load_sample= self.load_sample
        sparse = self.sparse
        format_list= self.format_list
        pool = mp.Pool(cpus)
        
        
        import time
        if vcf_file[-3:] == ".gz" or vcf_file[-4:] == ".bgz":
            infile = gzip.open(vcf_file, "rb")
            is_gzip = True
        else:
            infile = open(vcf_file, "r")
            is_gzip = False
        
        FixedINFO = {}
        contig_lines = []
        comment_lines = []
        var_ids, obs_ids, obs_dat = [], [], []
        count=0 #57077    
        for line in infile:
            count+=1
            # if count>10000:
            #     break
            if is_gzip:
                line = line.decode('utf-8')
            if line.startswith("#"):
                if line.startswith("##contig="):
                    contig_lines.append(line.rstrip())
                if line.startswith("#CHROM"):
                    if load_sample:
                        obs_ids = line.rstrip().split("\t")[9:]
                        for ob_id in obs_ids:
                            self.exclusive_donor_variants[ob_id]=set()
                    key_ids = line[1:].rstrip().split("\t")[:8]
                    for _key in key_ids:
                        FixedINFO[_key] = []
                else:
                    comment_lines.append(line.rstrip())
            else:
                pool.apply_async(self.load_sample_mp, args=([line,obs_ids,count,format_list]),callback=self.append_results)
                del line
        self.last_count=count
        pool.close()
        pool.join()
        
        output = self.combine_written_files()
        
        return output


"""Run CLI."""
parser = argparse.ArgumentParser(
    description="""
        Merges all the Celltype information in one file for Azimuth and Celltypist.
        """
)

parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s {version}'.format(version=__version__)
)

parser.add_argument(
    '-cpus', '--cpus',
    action='store',
    dest='cpus',
    required=True,
    help=''
)

parser.add_argument(
    '-cell_vcf', '--cell_vcf',
    action='store',
    dest='cell_vcf',
    required=True,
    help=''
)
# cell_assignments
parser.add_argument(
    '-cell_assignments', '--cell_assignments',
    action='store',
    dest='cell_assignments',
    required=True,
    help=''
)

parser.add_argument(
    '-donor_assignments', '--donor_assignments',
    action='store',
    dest='donor_assignments',
    required=True,
    help=''
)

parser.add_argument(
    '-gt_match_vcf', '--gt_match_vcf',
    action='store',
    dest='gt_match_vcf',
    required=True,
    help=''
)

parser.add_argument(
    '-expected_vcf', '--expected_vcf',
    action='store',
    dest='expected_vcf',
    required=True,
    help=''
)

cpus=10
cell_vcf = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/Cardinal_46291_Nov_29_2022/results_rsync2_copy/results/cellsnp/cellsnp_CRD_CMB13259712/cellSNP.cells.vcf.gz'
donor_assignments = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/Cardinal_46291_Nov_29_2022/results_rsync2_copy/results/gtmatch/CRD_CMB13259712/PiHAT_Stats_File_CRD_CMB13259712.csv'
gt_match_vcf = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/Cardinal_46291_Nov_29_2022/results_rsync2/results/subset_genotypes/Genotype_CRD_CMB13259712/InferedMerge_InferedGTMatched_CRD_CMB13259712.vcf.gz'
expected_vcf = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/Cardinal_46291_Nov_29_2022/results_rsync2/results/subset_genotypes/Genotype_CRD_CMB13259712/InferedMerge_InferedExpected_CRD_CMB13259712.vcf.gz'
cell_assignments = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/Cardinal_46291_Nov_29_2022/results_rsync2_copy/results/deconvolution/vireo_gt_fix/CRD_CMB13259712/GT_replace_donor_ids_false.tsv'


options = parser.parse_args()
cpus=int(options.cpus)
cell_vcf=options.cell_vcf
donor_assignments=options.donor_assignments
gt_match_vcf=options.gt_match_vcf
expected_vcf=options.expected_vcf
cell_assignments=options.cell_assignments

exclusive_donor_variants = {} #This is where results are populated when mp process i used.
curently_pushing =[] #this is a lock value to check if rhe curent field is updated so to avaid the race for update
All_Results={}
cell_concordance_table = {}

def append_results(result):
    # exclusive_donor_variants
    obs_with_gt= result[0]
    list_val_with_gt= result[1]
    idx = result[2]
    list_val = result[3]
    count = result[4]
    count11=0
    r = random.random()
    for ob_id in obs_with_gt:
        donor_loc_in_list = count11
        alleles = list_val_with_gt[donor_loc_in_list].split(':')[idx]
        if alleles!='.':
            ids = "_".join([list_val[x] for x in [0, 1, 3, 4]])
            donor_var = f"{ids}_{alleles}"
            while ob_id in curently_pushing:
               time.sleep(r*0.1)
            curently_pushing.append(ob_id)           
            exclusive_donor_variants[ob_id].add(donor_var)
            curently_pushing.remove(ob_id)
        count11+=1 
  
def find(lst, a):
    return [i for i, x in enumerate(lst) if x==a ]

def load_sample_mp(line,obs_ids,count,format_list):
    list_val = line.rstrip().split("\t") #[:5] #:8
    idx = find(list_val[8].split(':'),format_list[0])[0]
    if len(list_val[3]) > 1 or len(list_val[4]) > 1:
        # CURRENTLY DEALS ONLY WITH BIALELIC
        print(f'{idx} var not bialelic')
    else:
        list_val2 = list_val[9:]
        obs = pd.DataFrame(obs_ids)
        lv = pd.DataFrame(list_val2)
        lv_proc =lv[0].str.split(':').str[idx]
        gt_exists = lv_proc[lv_proc != '.']
        idx2 = gt_exists.index
        obs_with_gt = obs.loc[idx2.values]
        obs_with_gt = list(obs_with_gt[0].values)
        list_val_with_gt = lv.loc[idx2.values]
        list_val_with_gt = list(list_val_with_gt[0].values)
        random.seed(count)
        c = list(zip(obs_with_gt, list_val_with_gt))
        random.shuffle(c)
        obs_with_gt, list_val_with_gt = zip(*c)

    return [obs_with_gt,list_val_with_gt,idx,list_val,count]

def load_VCF_batch_paralel(vcf_file, biallelic_only=False, load_sample=True, sparse=True,
             format_list=None):
    """
    Load whole VCF file by utilising multiple cores to speed up loading of large cell files
    -------------------
    Initially designed to load VCF from cellSNP output, requiring 
    1) all variants have the same format list;
    2) a line starting with "#CHROM", with sample ids.
    If these two requirements are satisfied, this function also supports general
    VCF files, e.g., genotype for multiple samples.

    Note, it may take a large memory, please filter the VCF with bcftools first.
    """
    pool = mp.Pool(cpus)
    import time
    if vcf_file[-3:] == ".gz" or vcf_file[-4:] == ".bgz":
        infile = gzip.open(vcf_file, "rb")
        is_gzip = True
    else:
        infile = open(vcf_file, "r")
        is_gzip = False
    
    FixedINFO = {}
    contig_lines = []
    comment_lines = []
    var_ids, obs_ids, obs_dat = [], [], []
    
    # exclusive_donor_variants = {}
    count=0 #57077    
    for line in infile:
        count+=1
        # if count>10000:
        #     break
        if is_gzip:
            line = line.decode('utf-8')
        if line.startswith("#"):
            if line.startswith("##contig="):
                contig_lines.append(line.rstrip())
            if line.startswith("#CHROM"):
                if load_sample:
                    obs_ids = line.rstrip().split("\t")[9:]
                    for ob_id in obs_ids:
                        exclusive_donor_variants[ob_id]=set()
                key_ids = line[1:].rstrip().split("\t")[:8]
                for _key in key_ids:
                    FixedINFO[_key] = []
            else:
                comment_lines.append(line.rstrip())
        else:
            # #These lines activate will ignore the multiprocessing - was used to optimise the performance and check if there are no contestant for recording in directory.
            # line_count = load_sample_mp(line,obs_ids,count,format_list)       
            # append_results(line_count)
            ## Apply the multiprocessing strategy to load and process multiple lines simultaneously. 
            pool.apply_async(load_sample_mp, args=([line,obs_ids,count,format_list]),callback=append_results)
    pool.close()
    pool.join()
    return exclusive_donor_variants

def load_VCF_batch(vcf_file, biallelic_only=False, load_sample=True, sparse=True,
             format_list=None):
    """
    Load whole VCF file 
    -------------------
    Initially designed to load VCF from cellSNP output, requiring 
    1) all variants have the same format list;
    2) a line starting with "#CHROM", with sample ids.
    If these two requirements are satisfied, this function also supports general
    VCF files, e.g., genotype for multiple samples.

    Note, it may take a large memory, please filter the VCF with bcftools first.
    """
    
    if vcf_file[-3:] == ".gz" or vcf_file[-4:] == ".bgz":
        infile = gzip.open(vcf_file, "rb")
        is_gzip = True
    else:
        infile = open(vcf_file, "r")
        is_gzip = False
    
    FixedINFO = {}
    contig_lines = []
    comment_lines = []
    var_ids, obs_ids, obs_dat = [], [], []
    def find(lst, a):
        return [i for i, x in enumerate(lst) if x==a ]
    exclusive_donor_variants = {}
    count=0 #57077
    for line in infile:
        count+=1
        if is_gzip:
            line = line.decode('utf-8')
        if line.startswith("#"):
            if line.startswith("##contig="):
                contig_lines.append(line.rstrip())
            if line.startswith("#CHROM"):
                if load_sample:
                    obs_ids = line.rstrip().split("\t")[9:]
                    for ob_id in obs_ids:
                        exclusive_donor_variants[ob_id]=[]
                key_ids = line[1:].rstrip().split("\t")[:8]
                for _key in key_ids:
                    FixedINFO[_key] = []
            else:
                comment_lines.append(line.rstrip())
        else:
            list_val = line.rstrip().split("\t") #[:5] #:8
            idx = find(list_val[8].split(':'),format_list[0])[0]
            if biallelic_only:
                if len(list_val[3]) > 1 or len(list_val[4]) > 1:
                    continue
            if load_sample:
                # obs_dat.append(list_val[9:])
                list_val2 = list_val[9:]
                # len(list_val2)
                count11=0
                obs = pd.DataFrame(obs_ids)
                lv = pd.DataFrame(list_val2)
                lv_proc =lv[0].str.split(':').str[idx]
                gt_exists = lv_proc[lv_proc != '.']
                idx2 = gt_exists.index
                obs_with_gt = obs.loc[idx2.values]
                obs_with_gt = list(obs_with_gt[0].values)
                list_val_with_gt = lv.loc[idx2.values]
                list_val_with_gt = list(list_val_with_gt[0].values)
                for ob_id in obs_with_gt:
                    donor_loc_in_list = count11
                    alleles = list_val_with_gt[donor_loc_in_list].split(':')[idx]

                    if alleles!='.':
                        ids = "_".join([list_val[x] for x in [0, 1, 3, 4]])
                        donor_var = f"{ids}_{alleles}"
                        exclusive_donor_variants[ob_id].append(donor_var) 
                    count11+=1                      
    return exclusive_donor_variants

def load_VCF(vcf_file, biallelic_only=False, load_sample=True, sparse=True,
             format_list=None):
    """
    Load whole VCF file 
    -------------------
    Initially designed to load VCF from cellSNP output, requiring 
    1) all variants have the same format list;
    2) a line starting with "#CHROM", with sample ids.
    If these two requirements are satisfied, this function also supports general
    VCF files, e.g., genotype for multiple samples.

    Note, it may take a large memory, please filter the VCF with bcftools first.
    """
    if vcf_file[-3:] == ".gz" or vcf_file[-4:] == ".bgz":
        infile = gzip.open(vcf_file, "rb")
        is_gzip = True
    else:
        infile = open(vcf_file, "r")
        is_gzip = False
    
    FixedINFO = {}
    contig_lines = []
    comment_lines = []
    var_ids, obs_ids, obs_dat = [], [], []
    def find(lst, a):
        return [i for i, x in enumerate(lst) if x==a ]
    exclusive_donor_variants = {}
    
    for line in infile:
        if is_gzip:
            line = line.decode('utf-8')
        if line.startswith("#"):
            if line.startswith("##contig="):
                contig_lines.append(line.rstrip())
            if line.startswith("#CHROM"):
                if load_sample:
                    obs_ids = line.rstrip().split("\t")[9:]
                    for ob_id in obs_ids:
                        exclusive_donor_variants[ob_id]=[]
                key_ids = line[1:].rstrip().split("\t")[:8]
                for _key in key_ids:
                    FixedINFO[_key] = []
            else:
                comment_lines.append(line.rstrip())
        else:
            list_val = line.rstrip().split("\t") #[:5] #:8
            idx = find(list_val[8].split(':'),format_list[0])[0]
            if biallelic_only:
                if len(list_val[3]) > 1 or len(list_val[4]) > 1:
                    continue
            if load_sample:
                obs_dat.append(list_val[8:])
                count=0
                for ob_id in obs_ids:
                    donor_loc_in_list = count+9
                    ids = "_".join([list_val[x] for x in [0, 1, 3, 4]])
                    alleles = list_val[donor_loc_in_list].split(':')[idx]
                    donor_var = f"{ids}_{alleles}"
                    exclusive_donor_variants[ob_id].append(donor_var)
                    count+=1
                
            for i in range(len(key_ids)):
                FixedINFO[key_ids[i]].append(list_val[i])
            var_ids.append("_".join([list_val[x] for x in [0, 1, 3, 4]]))
    infile.close()

    RV = {}
    RV["variants"]  = var_ids
    RV["FixedINFO"] = FixedINFO
    RV["contigs"]   = contig_lines
    RV["comments"]  = comment_lines
    if load_sample:
        RV["samples"]   = obs_ids
        RV["GenoINFO"], RV["n_SNP_tagged"]  = parse_sample_info(
            obs_dat, sparse, format_list)
    return RV, exclusive_donor_variants

def parse_sample_info(sample_dat, sparse=True, format_list=None):
    """
    Parse genotype information for each sample
    Note, it requires the format for each variants to 
    be the same.
    """
    if sample_dat == [] or sample_dat is None:
        return None

    # require the same format for all variants
    format_all = [x[0].split(":") for x in sample_dat]
    if format_list is None:
        format_list = format_all[0]

    RV = {}
    n_SNP_tagged = np.zeros(len(format_list), np.int64)
    for _key in format_list:
        RV[_key] = []
    if sparse:
        ## sparse matrix requires all keys
        format_set_all = [set(x) for x in format_all]
        if format_set_all.count(set(format_list)) != len(format_all):
            print("Error: require the same format for all variants.")
            exit()

        RV['indices'] = []
        RV['indptr'] = [0]
        RV['shape'] = (len(sample_dat[0][1:]), len(sample_dat))
        missing_val = ":".join(["."] * len(format_list))
        
        cnt = 0
        for j in range(len(sample_dat)): #variant j
            _line = sample_dat[j]
            key_idx = [format_all[j].index(_key) for _key in format_list]
            for i in range(len(_line[1:])): #cell i
                if _line[i+1] == missing_val or _line[i+1] == ".":
                    continue
                _line_key = _line[i+1].split(":")
                for k in range(len(format_list)):
                    RV[format_list[k]].append(_line_key[key_idx[k]])

                cnt += 1
                RV['indices'].append(i)
                n_SNP_tagged += 1
            RV['indptr'].append(cnt)
    else:
        for j in range(len(sample_dat)): #variant j
            _line = sample_dat[j]
            _line_split = [x.split(":") for x in _line[1:]]
            for il, _key in enumerate(format_list):
                if _key in format_all[j]:
                    k = format_all[j].index(_key)
                    _line_key = [x[k] for x in _line_split]
                    RV[_key].append(_line_key)
                    n_SNP_tagged[il] += 1
                else:
                    RV[_key].append(["."] * len(_line_split))

    # Check if format tags are well convered
    idx_low_tag = np.where(n_SNP_tagged < (0.1 * len(sample_dat)))[0]
    if len(idx_low_tag) > 0:
        print('[vireo] Warning: too few variants with tags!',
              '\t'.join([format_list[k] + ": " + str(n_SNP_tagged[k])
                         for k in range(len(format_list))]))
    
    return RV, n_SNP_tagged

def concordance_table(donor_cell_vcf,donor_genotype_replicated,donor_in_question,donor_gt_match):
    # Now that we have correctly named 2 gtfiles we can calculate the concordances/discordances
    # Load the Genotypes in vcf format as per vireo  and calculate the concordances as per Hail.
    biallelic_dataset_donor_cell_vcf = donor_cell_vcf.filter_rows(hl.len(donor_cell_vcf.alleles) == 2)
    biallelic_dataset_donor_genotype_replicated = donor_genotype_replicated.filter_rows(hl.len(donor_genotype_replicated.alleles) == 2)
    samples =  hl.concordance(biallelic_dataset_donor_cell_vcf, biallelic_dataset_donor_genotype_replicated)
    samples2 = samples.to_pandas(flatten=False)
    dataset= []
    for i in range(0,len(samples2)):
        id_2=samples2.loc[i,'s']
        total_concordant = samples2.iloc[i]['concordance'][2][2]+samples2.iloc[i]['concordance'][3][3]+samples2.iloc[i]['concordance'][4][4]
        total_discordant = sum([sum(s[2:]) for s in samples2.iloc[i]['concordance'][2:]]) - total_concordant
        percent_discordant = total_discordant/(total_concordant+total_discordant)*100
        dataset.append({'donor_in_question':donor_in_question,'donor_gt_match':donor_gt_match,'sample':id_2,'total_concordant':total_concordant,'total_discordant':total_discordant,'percent_discordant':percent_discordant})

    Data2 = pd.DataFrame(dataset) 
    del biallelic_dataset_donor_cell_vcf
    del biallelic_dataset_donor_genotype_replicated
    return Data2

def norm_genotypes(expected_vars):
    expected_vars = pd.DataFrame(expected_vars)
    expected_vars['ids'] = expected_vars[0].str.split("_").str[:-1].str.join('_')
    expected_vars['vars'] = expected_vars[0].str.split("_").str[-1].str.join('')
    expected_vars['vars'] = expected_vars['vars'].str.replace('|','/',regex=False)
    expected_vars = expected_vars[expected_vars['vars']!='./.']
    expected_vars.loc[expected_vars['vars']=='0/1','vars']='1/0'
    expected_vars['combo']= expected_vars['ids']+'_'+expected_vars['vars']
    return expected_vars

def retrieve_concordant_discordant_sites(expected_vars,cell_vars):
    # This function has been inspired by Hails Concordance implementations, however hail has a pitfall that it performs a lot of other stuff under hood and requires intermediate sorting operations.
    # Since the single cell calculations requires concordance calculations per cell this becomes very computationally heavyon Hail, hence we have implemented concordance calculations here as part of the pipeline.
    # Author: M.Ozols
    expected_vars_norm = norm_genotypes(expected_vars)
    cell_vars_norm = norm_genotypes(cell_vars)
    Total_Overlappin_sites = set(expected_vars_norm['ids']).intersection(set(cell_vars_norm['ids']))
    expected_vars2 = expected_vars_norm[expected_vars_norm['ids'].isin(Total_Overlappin_sites)]
    cell_vars2 = cell_vars_norm[cell_vars_norm['ids'].isin(Total_Overlappin_sites)]
    Concordant_Sites = set(cell_vars2['combo']).intersection(set(expected_vars2['combo']))
    Discodrant_sites = set(cell_vars2['combo'])-set(expected_vars2['combo'])
    
    return Concordant_Sites, Discodrant_sites, Total_Overlappin_sites

def donor_exclusive_sites(exclusive_don_variants2):
    # Here we generate a function for determining the sites that are donor exclusive
    donor_distinct_sites = {}
    for col1 in exclusive_don_variants2.keys():
        comparisons =[]
        to_compare = []
        for col2 in exclusive_don_variants2.keys():
            if col1==col2:
                # we set this as the unique entry
                # print('1')
                to_compare = set(exclusive_don_variants2[col2])
            else:
                # We combine all the variants in one list
                comparisons+=list(exclusive_don_variants2[col2])
                # print('2')
        # print('comparison')  
        comparisons_all = set(comparisons) 
        comparisons_all_norm = norm_genotypes(comparisons_all)
        comparisons_all=set(comparisons_all_norm['combo'])
        
        to_compare = set(to_compare)
        to_compare_norm = norm_genotypes(to_compare)
        to_compare=set(to_compare_norm['combo'])
        # Make sure we account for hap types - phased/unphased 
        distinct_donor_sites = to_compare - comparisons_all
        donor_distinct_sites[col1]=distinct_donor_sites
        # Perform the distinct set function.
    return donor_distinct_sites

def calculate_concordances(all_cells_to_perform_operation_on,donor_cell_vcf,donor_genotype_replicated,donor_in_question,donor_gt_match,count):
    Cells_to_keep =hl.literal(set(all_cells_to_perform_operation_on))
    donor_single_cell_vcf = donor_cell_vcf.filter_cols(Cells_to_keep.contains(donor_cell_vcf['s']))
    del Cells_to_keep
    # donor_genotype_replicated = donor_genotype_replicated.checkpoint(f'/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/cellsnp_optimisations/checkpoints/donor_checkpoint_{count}', overwrite=True)
    D2 = concordance_table(donor_single_cell_vcf,donor_genotype_replicated,donor_in_question,donor_gt_match)
    del donor_single_cell_vcf
    D2.to_csv(f'/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/cellsnp_optimisations/{count}_out.csv',sep='\t')
    # Data_All= pd.concat([Data_All,D2])
    return {'idx':count,'Result':D2}

    # if we can load n cells at a time from the cellvcf file then we can see if the concordant sites are correctly picked up and which sites are concordant and which are not.

def append_results_cell_concordances(result):
    # print(result)
    try:
        percent_concordant = result[2]/(result[3]+result[2])*100
    except:
        percent_concordant = 0
    
    try:
        percent_discordant = result[3]/(result[3]+result[2])*100
    except:
        percent_discordant = 0
    
    cell_concordance_table[f'{result[0]} --- {result[1]}'] = {'GT 1':result[0],
                                                            'GT 2':result[1],
                                                            'Nr_Concordant':result[2],
                                                            'Nr_Discordant':result[3],
                                                            'Percent Concordant':percent_concordant,
                                                            'Percent Discordant':percent_discordant,
                                                            'NrTotal_Overlapping_sites between two gewnotypes':result[4],
                                                            'Nr_donor_distinct_sites_within_pool_individuals':result[6],
                                                            'Number_of_sites_that_are_donor_concordant_and_exclusive':result[5],
                                                            }
    

def concordance_dable_production(expected_vars,cell_vars,cell1,donor_gt_match,dds):
    Nr_donor_distinct_sites = len(dds)
    Concordant_Sites, Discodrant_sites, Total_Overlappin_sites = retrieve_concordant_discordant_sites(expected_vars,cell_vars)
    Nr_Concordant = len(Concordant_Sites)
    Nr_Discordant = len(Discodrant_sites)
    Nr_Total_Overlapping_sites = len(Total_Overlappin_sites)
    Number_of_sites_that_are_donor_concordant_and_exclusive = len(set(dds).intersection(set(Concordant_Sites)))
    return [cell1,donor_gt_match,Nr_Concordant,Nr_Discordant,Nr_Total_Overlapping_sites,Number_of_sites_that_are_donor_concordant_and_exclusive,Nr_donor_distinct_sites]

def conc_table(donor_assignments_table,cell_assignments_table,exclusive_don_variants,exclusive_cell_variants):
    pool = mp.Pool(cpus)
    for i,row1 in donor_assignments_table.iterrows():
        donor_in_question = row1['donor_query']
        donor_gt_match = row1['donor_gt']
        # print(donor_gt_match)
        if (donor_gt_match=='NONE'):
            continue
        Cells_to_keep_pre = list(set(cell_assignments_table.loc[cell_assignments_table['donor_id']==donor_in_question,'cell']))
        try:
            # Now we subset this down to each of the uniqie variants per donor and check which of the concordant sites are exclusive to donor.
            dds = donor_distinct_sites[donor_gt_match]
        except:
            continue
        
        count = 0
        for cell1 in Cells_to_keep_pre:
            count+=1
            # if count>100:
            #     break
            # print(count)
            expected_vars = exclusive_don_variants[donor_gt_match]
            cell_vars = exclusive_cell_variants[cell1]
            cell_concordance_table[f'{cell1} --- {donor_gt_match}']={}
            pool.apply_async(concordance_dable_production, args=([expected_vars,cell_vars,cell1,donor_gt_match,dds]),callback=append_results_cell_concordances)          
    pool.close()
    pool.join()
    return cell_concordance_table

if __name__ == "__main__":
 
    print('---Genotype loader init----')    
    loader2 = VCF_Loader(gt_match_vcf, biallelic_only=True,
                    sparse=False, format_list=['GT'])
    GT_Matched_variants = loader2.load_VCF_batch_paralel()
    del loader2
    
     
    print('---Lets load cell vcf----')
    tic = time.perf_counter()
    loader1 = VCF_Loader(cell_vcf, biallelic_only=True,
                        sparse=False, format_list=['GT'])
    exclusive_cell_variants = loader1.load_VCF_batch_paralel()
    del loader1
    toc = time.perf_counter()
    print(f"Loadiong took {toc - tic:0.4f} seconds")
    # exclusive_cell_variants = load_VCF_batch_paralel(cell_vcf, biallelic_only=True,
    #                     sparse=False, format_list=['GT'])
    
    print('---Cell VCF file loaded----')
    donor_assignments_table = pd.read_csv(donor_assignments)
    cell_assignments_table = pd.read_csv(cell_assignments,sep='\t')
    

    print('---Variant1 files loaded----')
    loader3 = VCF_Loader(expected_vcf, biallelic_only=True,
                    sparse=False, format_list=['GT'])
    GT_Expected_variants = loader3.load_VCF_batch_paralel()
    del loader3
    
    
    # GT_Matched_variants2 = load_VCF_batch(gt_match_vcf, biallelic_only=True,
    #                         sparse=False, format_list=['GT'])
    # GT_Expected_variants = load_VCF_batch(expected_vcf, biallelic_only=True,
    #                         sparse=False, format_list=['GT'])
    print('---Variant2 files loaded----')
    exclusive_don_variants = GT_Expected_variants.keys()
    content = [x for x in exclusive_don_variants if not x.startswith('donor')]
    GT_Expected_variants = {key: GT_Expected_variants[key] for key in content}
    
    exclusive_don_variants = GT_Matched_variants.keys()
    content = [x for x in exclusive_don_variants if not x.startswith('donor')]
    GT_Matched_variants = {key: GT_Matched_variants[key] for key in content}
    
    exclusive_don_variants = GT_Expected_variants
    for key in GT_Matched_variants.keys():
        if key in exclusive_don_variants.keys():
            _=''
        else:
            exclusive_don_variants[key]=GT_Matched_variants[key]
        
    donor_distinct_sites = donor_exclusive_sites(exclusive_don_variants)
    print('---donor_distinct_sites calculated----')
    cell_concordance_table = conc_table(donor_assignments_table,cell_assignments_table,exclusive_don_variants,exclusive_cell_variants)
    result = pd.DataFrame(cell_concordance_table).T
    result.to_csv('cell_concordance_table.tsv',sep='\t')
    print('Processing Done')
    
