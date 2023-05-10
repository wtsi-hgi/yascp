#!/usr/bin/env python3

__date__ = '2023-05-10'
__version__ = '0.0.1'
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


class Concordances:
        def __init__(self, donor_assignments_table,cell_assignments_table,exclusive_don_variants,exclusive_cell_variants,exclusive_cell_variants_dp,donor_distinct_sites):
            self.reset()
            self.donor_assignments_table=donor_assignments_table
            self.cell_assignments_table=cell_assignments_table
            self.exclusive_don_variants=exclusive_don_variants
            self.exclusive_cell_variants=exclusive_cell_variants
            self.exclusive_cell_variants_dp=exclusive_cell_variants_dp
            self.donor_distinct_sites=donor_distinct_sites
            self.record_dict={}

        def norm_genotypes(self,expected_vars):
            expected_vars = pd.DataFrame(expected_vars)
            if len(expected_vars) > 0:
                split_str=expected_vars[0].str.split("_")
                expected_vars['ids'] = split_str.str[0]+'_'+split_str.str[1]+'_'+split_str.str[2]+'_'+split_str.str[3]
                expected_vars['pos'] = split_str.str[0]+'_'+split_str.str[1]
                expected_vars['vars'] = split_str.str[-1]
                expected_vars['vars'] = expected_vars['vars'].str.replace('|','/',regex=False)
                expected_vars = expected_vars[expected_vars['vars']!='./.']
                expected_vars.loc[expected_vars['vars']=='0/1','vars']='1/0'
                expected_vars['combo']= expected_vars['ids']+'_'+expected_vars['vars']
            return expected_vars
        
        def reset(self):
            self.cell_concordance_table ={}


        def get_strict_discordance(self, snp_gtypes, cellsnp_gtypes):
            '''
            take a list of SNP array genotypes and a list of cellSNP genotypes, return counts of truly discordant 
            sites and relaxed concordant sites
            1) If you have 1/1 om SNP array you can not get a 0/1 or 0/0 genotype
            2) if you have a 0/0 you can not get a 1/1 or 0/1
            3) if you genotype is 0/1 you can get all copies: 0/0 . 0/1. 1/1
            So - each obversed cellsnp allele must be in the array SNP gtype
            '''
            true_discordant = 0
            relaxed_concordant = 0

            for i in range(0, len(snp_gtypes)):
                discordant = False
                snp_alleles = [snp_gtypes[i][-3], snp_gtypes[i][-1]]
                cellsnp_alleles = [cellsnp_gtypes[i][-3], cellsnp_gtypes[i][-1]]
                snp_alleles_set = set(snp_alleles)
                cellsnp_alleles_set = set(cellsnp_alleles)
                
                snp_var = snp_gtypes[i][:-3]
                cellsnp_var = cellsnp_gtypes[i][:-3]

                if not cellsnp_var == snp_var:
                    print("Error with strict discordance calculations: " + snp_gtypes[i] + " " + cellsnp_gtypes[i])
                    exit(1)
                else:
                    for allele in cellsnp_alleles_set:
                        if not allele in snp_alleles_set:#if a cellSNP allele is found that is not in the array data this is discordant
                            discordant = True
                
                if discordant == True:
                    true_discordant+=1
                else:
                    relaxed_concordant+=1

            return true_discordant, relaxed_concordant
        

        def retrieve_concordant_discordant_sites(self,expected_vars_norm,cell_vars):
            # This function has been inspired by Hails Concordance implementations, however hail has a pitfall that it performs a lot of other stuff under hood and requires intermediate sorting operations.
            # Since the single cell calculations requires concordance calculations per cell this becomes very computationally heavy on Hail, hence we have implemented concordance calculations here as part of the pipeline.
            # Author: M.Ozols
            
            cell_vars_norm = self.norm_genotypes(cell_vars)
            if len(cell_vars_norm) > 0:
                Total_Overlappin_sites = set(expected_vars_norm['ids']).intersection(set(cell_vars_norm['ids']))
                expected_vars2 = expected_vars_norm[expected_vars_norm['ids'].isin(Total_Overlappin_sites)]
                cell_vars2 = cell_vars_norm[cell_vars_norm['ids'].isin(Total_Overlappin_sites)]
                Concordant_Sites = set(cell_vars2['combo']).intersection(set(expected_vars2['combo']))
                Discodrant_sites = set(cell_vars2['combo'])-set(expected_vars2['combo'])
                disc = pd.DataFrame(Discodrant_sites,columns=['combo_x'])
                df_cd = pd.merge(cell_vars2, expected_vars2, how='inner', on = 'pos')
                disc2= pd.merge(disc, df_cd, how='inner', on = 'combo_x')
                disc2['expected_retrieved'] = disc2['0_x']+'::'+disc2['0_y']
                disc_sites = ';'.join(disc2['expected_retrieved'])
                #find truly discordant sites
                true_discordant_count, relaxed_concordant_count = self.get_strict_discordance(disc2['0_y'], disc2['0_x'])
            else:
                Total_Overlappin_sites = set()
                Concordant_Sites = set()
                Discodrant_sites = set()
                disc_sites = ''
                true_discordant_count = 0
                relaxed_concordant_count = 0

            return Concordant_Sites, Discodrant_sites, Total_Overlappin_sites, disc_sites,cell_vars_norm, true_discordant_count, relaxed_concordant_count
        

        def set_results(self,to_set,id):
            # Recod to disk to save the loading mmeory time.
            with open(f'tmp_{id}.pkl', 'wb') as f:
                pickle.dump(to_set, f)
            self.record_dict[id]=f'tmp_{id}.pkl'
        
        def append_results_cell_concordances(self,result):
            count=result[9]
            try:
                percent_concordant = result[2]/(result[3]+result[2])*100
            except:
                percent_concordant = 0
            
            try:
                percent_discordant = result[3]/(result[3]+result[2])*100
            except:
                percent_discordant = 0

            try:
                percent_relaxed_concordant = result[4]/(result[4]+result[5])*100
            except:
                percent_relaxed_concordant = 0
            
            try:
                percent_strict_discordant = result[5]/(result[4]+result[5])*100
            except:
                percent_strict_discordant = 0

            try:
                percent_concordant_dp = result[11]/(result[12]+result[11])*100
            except:
                percent_concordant_dp = 0

            try:
                percent_discordant_dp = result[12]/(result[12]+result[11])*100
            except:
                percent_discordant_dp = 0

            try:
                percent_relaxed_concordant_dp = result[13]/(result[13]+result[14])*100
            except:
                percent_relaxed_concordant_dp = 0
            
            try:
                percent_strict_discordant_dp = result[14]/(result[13]+result[14])*100
            except:
                percent_strict_discordant_dp = 0


                # self.reset_c()   
            print(count)
            self.cell_concordance_table[f'{result[0]} --- {result[1]}'] = {'GT 1':result[0],
                                                                    'GT 2':result[1],
                                                                    'Nr_Concordant':result[2],
                                                                    'Nr_Discordant':result[3],
                                                                    'Nr_Relaxed_concordant':result[4],
                                                                    'Nr_strict_discordant':result[5],
                                                                    'Percent Concordant':percent_concordant,
                                                                    'Percent Discordant':percent_discordant,
                                                                    'Percent_relaxed_concordant': percent_relaxed_concordant,
                                                                    'Percent_strict_discordant': percent_strict_discordant,
                                                                    'NrTotal_Overlapping_sites_between_two_genotypes':result[6],
                                                                    'Nr_donor_distinct_sites_within_pool_individuals':result[8],
                                                                    'Number_of_sites_that_are_donor_concordant_and_exclusive':result[7],
                                                                    'Discordant_Site_Identities':result[10],
                                                                    'Nr_Concordant_dp':result[11],
                                                                    'Nr_Discordant_dp':result[12],
                                                                    'Nr_relaxed_concordant_dp':result[13],
                                                                    'Nr_strict_discordant_dp':result[14],
                                                                    'Percent_Concordant_dp':percent_concordant_dp,
                                                                    'Percent_Discordant_dp':percent_discordant_dp,
                                                                    'Percent_relaxed_concordant_dp': percent_relaxed_concordant_dp,
                                                                    'Percent_strict_discordant_dp': percent_strict_discordant_dp,
                                                                    'NrTotal_Overlapping_sites_between_two_genotypes_dp':result[15],
                                                                    'Number_of_sites_that_are_donor_concordant_and_exclusive_dp':result[16],
                                                                    'Discordant_Site_Identities_dp':result[17]
                                                                    }   
            
            if (count % 200 == 0):
                print(f'recording and resetting memory {count}')
                # self.record_dict[count]=self.exclusive_donor_variants
                self.set_results(self.cell_concordance_table,count)
                self.reset()  
            _=""
        
        def combine_written_files(self):
            to_export = self.cell_concordance_table
            for val1 in self.record_dict.values():
                # here remove the int files.
                print(f"merging temp file: {val1}")
                with open(val1, 'rb') as f:
                    loaded_dict = pickle.load(f)
                    for k1 in loaded_dict.keys():
                        to_export[k1]=loaded_dict[k1]
                os.remove(val1)
            return to_export
        
        
        def conc_table(self):
            donor_assignments_table=self.donor_assignments_table
            cell_assignments_table=self.cell_assignments_table
            exclusive_don_variants=self.exclusive_don_variants
            exclusive_cell_variants= self.exclusive_cell_variants
            exclusive_cell_variants_dp = self.exclusive_cell_variants_dp
            
            pool = mp.Pool(cpus)
            count = 0
            for i,row1 in donor_assignments_table.iterrows():
                donor_in_question = row1['donor_query']
                donor_gt_match = row1['donor_gt']
                # print(donor_gt_match)
                if (donor_gt_match=='NONE'):
                    continue
                Cells_to_keep_pre = list(set(cell_assignments_table.loc[cell_assignments_table['donor_id']==donor_in_question,'cell']))
                expected_vars = exclusive_don_variants[donor_gt_match]
                expected_vars_norm = self.norm_genotypes(expected_vars)
                try:
                    # Now we subset this down to each of the uniqie variants per donor and check which of the concordant sites are exclusive to donor.
                    dds = self.donor_distinct_sites[donor_gt_match]
                except:
                    continue
                
                for cell1 in Cells_to_keep_pre:
                    count+=1
                    # if count>800:
                    #     break
                    cell_vars = exclusive_cell_variants[cell1]
                    cell_vars_dp = exclusive_cell_variants_dp[cell1]

                    self.cell_concordance_table[f'{cell1} --- {donor_gt_match}']={}
                    # pool.apply_async(self.concordance_dable_production, args=([expected_vars_norm,cell_vars,cell1,donor_gt_match,dds,count]),callback=self.append_results_cell_concordances)          
                    result1 = self.concordance_dable_production(expected_vars_norm,cell_vars,cell_vars_dp,cell1,donor_gt_match,dds,count)
                    self.append_results_cell_concordances(result1)
                    
            pool.close()
            pool.join()
            output = self.combine_written_files()
            return output
        
        def concordance_dable_production(self,expected_vars_norm,cell_vars,cell_vars_dp,cell1,donor_gt_match,dds,count):
            Nr_donor_distinct_sites = len(dds)
            Concordant_Sites, Discodrant_sites, Total_Overlappin_sites,discordant_sites,cell_vars_norm, Nr_strict_discordant, relaxed_concordant_count = self.retrieve_concordant_discordant_sites(expected_vars_norm,cell_vars)
            Nr_Concordant = len(Concordant_Sites)
            Nr_Relaxed_concordant = Nr_Concordant + relaxed_concordant_count
            Nr_Discordant = len(Discodrant_sites)
            Nr_Total_Overlapping_sites = len(Total_Overlappin_sites)
            Number_of_sites_that_are_donor_concordant_and_exclusive = len(set(dds).intersection(set(Concordant_Sites)))
            Number_of_sites_in_cellsnp_but_not_in_reference = set(cell_vars_norm['pos'])-set(expected_vars_norm['pos'])

            Concordant_Sites_dp, Discodrant_sites_dp, Total_Overlappin_sites_dp,discordant_sites_dp,cell_vars_norm_dp, Nr_strict_discordant_dp, relaxed_concordant_count_dp = self.retrieve_concordant_discordant_sites(expected_vars_norm,cell_vars_dp)
            Nr_Concordant_dp = len(Concordant_Sites_dp)
            Nr_Relaxed_concordant_dp = Nr_Concordant_dp + relaxed_concordant_count_dp
            Nr_Discordant_dp = len(Discodrant_sites_dp)
            Nr_Total_Overlapping_sites_dp = len(Total_Overlappin_sites_dp)
            Number_of_sites_that_are_donor_concordant_and_exclusive_dp = len(set(dds).intersection(set(Concordant_Sites_dp)))
            Number_of_sites_in_cellsnp_but_not_in_reference = set(cell_vars_norm['pos'])-set(expected_vars_norm['pos'])

            return [cell1,donor_gt_match,Nr_Concordant,Nr_Discordant,Nr_Relaxed_concordant, Nr_strict_discordant, Nr_Total_Overlapping_sites,
                    Number_of_sites_that_are_donor_concordant_and_exclusive, Nr_donor_distinct_sites,count,discordant_sites,
                    Nr_Concordant_dp,Nr_Discordant_dp, Nr_Relaxed_concordant_dp, Nr_strict_discordant_dp, Nr_Total_Overlapping_sites_dp,
                    Number_of_sites_that_are_donor_concordant_and_exclusive_dp,discordant_sites_dp]
        
        
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
        if (count % 200 == 0):
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

parser.add_argument(
    '-cell_vcf_dp', '--cell_vcf_dp',
    action='store',
    dest='cell_vcf_dp',
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
cell_vcf_dp=options.cell_vcf_dp#cell VCF data filtered by DP
donor_assignments=options.donor_assignments
gt_match_vcf=options.gt_match_vcf
expected_vcf=options.expected_vcf
cell_assignments=options.cell_assignments

exclusive_donor_variants = {} #This is where results are populated when mp process i used.
curently_pushing =[] #this is a lock value to check if rhe curent field is updated so to avaid the race for update
All_Results={}
cell_concordance_table = {}

  
def find(lst, a):
    return [i for i, x in enumerate(lst) if x==a ]


def norm_genotypes(expected_vars):
    expected_vars = pd.DataFrame(expected_vars)
    expected_vars['ids'] = expected_vars[0].str.split("_").str[:-1].str.join('_')
    expected_vars['pos'] = expected_vars[0].str.split("_").str[:2].str.join('_')
    expected_vars['vars'] = expected_vars[0].str.split("_").str[-1].str.join('')
    expected_vars['vars'] = expected_vars['vars'].str.replace('|','/',regex=False)
    expected_vars = expected_vars[expected_vars['vars']!='./.']
    expected_vars.loc[expected_vars['vars']=='0/1','vars']='1/0'
    expected_vars['combo']= expected_vars['ids']+'_'+expected_vars['vars']
    return expected_vars


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

debug=False

if __name__ == "__main__":
    donor_assignments_table = pd.read_csv(donor_assignments)
    cell_assignments_table = pd.read_csv(cell_assignments,sep='\t')
    
    if debug:
        with open('tmp_GT_Expected_variants.pkl', 'rb') as f:
            GT_Expected_variants = pickle.load(f)
        with open('tmp_GT_Matched_variants.pkl', 'rb') as f:
            GT_Matched_variants = pickle.load(f)   
        with open('tmp_exclusive_cell_variants.pkl', 'rb') as f:
            exclusive_cell_variants = pickle.load(f) 
        with open('tmp_donor_distinct_sites.pkl', 'rb') as f:
            donor_distinct_sites = pickle.load(f)  
        with open('tmp_exclusive_don_variants.pkl', 'rb') as f:
            exclusive_don_variants = pickle.load(f) 
        with open('tmp_exclusive_cell_variants_dp.pkl', 'rb') as f:
            exclusive_cell_variants_dp = pickle.load(f) 
    else:  
        print('---Genotype loader init----')    
        loader2 = VCF_Loader(gt_match_vcf, biallelic_only=True,
                        sparse=False, format_list=['GT'])
        GT_Matched_variants = loader2.load_VCF_batch_paralel()
        del loader2

        with open(f'tmp_GT_Matched_variants.pkl', 'wb') as f:
            pickle.dump(GT_Matched_variants, f)
        
        print('---Lets load cell vcf----')
        tic = time.perf_counter()
        loader1 = VCF_Loader(cell_vcf, biallelic_only=True,
                            sparse=False, format_list=['GT'])
        exclusive_cell_variants = loader1.load_VCF_batch_paralel()
        del loader1

        toc = time.perf_counter()
            
        with open(f'tmp_exclusive_cell_variants.pkl', 'wb') as f:
            pickle.dump(exclusive_cell_variants, f)
        print(f"Loadiong took {toc - tic:0.4f} seconds")

        # exclusive_cell_variants = load_VCF_batch_paralel(cell_vcf, biallelic_only=True,
        #                     sparse=False, format_list=['GT'])
        
        print('---Cell VCF file loaded----')
      
        print('---Lets load DP filtered cell vcf----')
        loader4 = VCF_Loader(cell_vcf_dp, biallelic_only=True,
                            sparse=False, format_list=['GT'])
        exclusive_cell_variants_dp = loader4.load_VCF_batch_paralel()
        del loader4

        with open(f'tmp_exclusive_cell_variants_dp.pkl', 'wb') as f:
            pickle.dump(exclusive_cell_variants_dp, f)

        print('---Variant1 files loaded----')
        loader3 = VCF_Loader(expected_vcf, biallelic_only=True,
                        sparse=False, format_list=['GT'])
        GT_Expected_variants = loader3.load_VCF_batch_paralel()
        del loader3

        with open(f'tmp_GT_Expected_variants.pkl', 'wb') as f:
            pickle.dump(GT_Expected_variants, f)
        
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
        
        with open(f'tmp_exclusive_don_variants.pkl', 'wb') as f:
            pickle.dump(exclusive_don_variants, f)
        donor_distinct_sites = donor_exclusive_sites(exclusive_don_variants)
        with open(f'tmp_donor_distinct_sites.pkl', 'wb') as f:
            pickle.dump(donor_distinct_sites, f)
        
    print('---donor_distinct_sites calculated----')
    
    conc1 = Concordances(donor_assignments_table,cell_assignments_table,exclusive_don_variants,exclusive_cell_variants,exclusive_cell_variants_dp,donor_distinct_sites)
    cell_concordance_table = conc1.conc_table()
    
    # cell_concordance_table = conc_table(donor_assignments_table,cell_assignments_table,exclusive_don_variants,exclusive_cell_variants)
    result = pd.DataFrame(cell_concordance_table).T
    result.to_csv('cell_concordance_table.tsv',sep='\t')
    print('Processing Done')
    
