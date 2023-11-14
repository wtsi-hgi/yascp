#!/usr/bin/env python3

#take cellSNP VCF and genotype VCF for the donors in a pool
# for each cell in the cellSNP VCF identify discordant sites (using the relaxed concordance)
# look for these sites in genotypes of all members of the pool
# output:
#   cell id
#   assigned donor
#   cohort of assigned donor
#   number of discordant sites
#   total AD over discordant sites
#   list of donors in the pool, how many of the discordant sites are found in the donor, cohort each belongs to
#   list of discordant sites

__date__ = '2023-14-11'
__version__ = '0.0.1'
import argparse
import sys
import importlib.util
import pickle 
import pandas as pd
import gzip
import random
import numpy as np
import time
import multiprocessing as mp
from multiprocessing import Lock
import logging
import os
import gzip
import time
pd.options.mode.chained_assignment = None

class Concordances:
    
    def reset(self):
        self.cell_concordance_table ={}

    def reset2(self):
        self.other_donor_comp =[] 
    
    def __init__(self, donor_assignments_table,cell_assignments_table,exclusive_don_variants,exclusive_cell_variants,donor_distinct_sites,informative_sites, uninformative_sites):
        self.reset()
        self.reset2()
        self.donor_assignments_table=donor_assignments_table
        self.cell_assignments_table=cell_assignments_table
        self.exclusive_don_variants=exclusive_don_variants
        self.exclusive_cell_variants=exclusive_cell_variants
        self.donor_distinct_sites=donor_distinct_sites
        self.informative_sites = informative_sites
        self.uninformative_sites = uninformative_sites
        self.record_dict={}

    def norm_genotypes(self,expected_vars):
        expected_vars = pd.DataFrame(expected_vars)
        if len(expected_vars) > 0:
            split_str=expected_vars[0].str.split("_")
            expected_vars['ids'] = split_str.str[0]+'_'+split_str.str[1]+'_'+split_str.str[2]+'_'+split_str.str[3]
            expected_vars['pos'] = split_str.str[0]+'_'+split_str.str[1]
            expected_vars['vars'] = split_str.str[4]
            expected_vars['vars'] = expected_vars['vars'].str.replace('|','/',regex=False)
            expected_vars = expected_vars[expected_vars['vars']!='./.']
            expected_vars.loc[expected_vars['vars']=='0/1','vars']='1/0'
            expected_vars['combo']= expected_vars['ids']+'_'+expected_vars['vars']
        return expected_vars



    def get_strict_discordance(self, expected_vars, cell_vars):
        '''
        take a list of SNP array genotypes and a list of cellSNP genotypes, return counts of truly discordant 
        sites and relaxed concordant sites and list of discordant sites1
        1) If you have 1/1 on SNP array you can not get a 0/1 or 0/0 genotype
        2) if you have a 0/0 you can not get a 1/1 or 0/1
        3) if you genotype is 0/1 you can get all copies: 0/0 . 0/1. 1/1
        So - each obversed cellsnp allele must be in the array SNP gtype
        '''
        snp_gtypes = expected_vars[0]
        cellsnp_gtypes = cell_vars[0]
        true_discordant = 0
        relaxed_concordant = 0
        relaxed_concordant_informative = 0
        relaxed_concordant_informative_ids = []
        relaxed_concordant_uninformative_ids = []
        true_discordant_uninformative_ids =[]
        true_discordant_informative_ids=[]
        relaxed_concordant_uninformative = 0
        true_discordant_informative = 0
        true_discordant_uninformative = 0
        discordant_vars = []
        concordant_vars = []
        subset_informative_concordant = 0
        subset_informative_discordant = 0

        #print(self.uninformative_sites)
        #print(self.informative_sites)

        #create sets of the ids (chrom, pos, ref, alt) in each set of genotypes. Filter to the ids present in both 
        #then filter to informative and uninformative. If uninformative >0 then create a subset of informative
        # with the same number of vars (at random)
        split_snp_gts=snp_gtypes.str.split("_")
        snp_gtypes_ids = set(split_snp_gts.str[0]+'_'+split_snp_gts.str[1]+'_'+split_snp_gts.str[2]+'_'+split_snp_gts.str[3])

        split_cellsnp_gts=cellsnp_gtypes.str.split("_")
        cellsnp_gtypes_ids = set(split_cellsnp_gts.str[0]+'_'+split_cellsnp_gts.str[1]+'_'+split_cellsnp_gts.str[2]+'_'+split_cellsnp_gts.str[3])

        shared_gts = snp_gtypes_ids.intersection(cellsnp_gtypes_ids)

        shared_informative = shared_gts.intersection(self.informative_sites)
        shared_uninformative = shared_gts.intersection(self.uninformative_sites)
        # print("shared informative " + str(len(shared_informative)))
        # print("shared uninformative " + str(len(shared_uninformative)))

        #store the numbers of informative and uninformative sites shared between cellSNP and gt data as these
        #are the sites used for concordance
        self.informative_covered = len(shared_informative)
        self.uninformative_covered = len(shared_uninformative)

        if len(shared_uninformative) > 0:
            #print(len(shared_uninformative))
            # print(len(shared_informative))
            if len(shared_uninformative) <= len(shared_informative):
                informative_subset = set(random.sample(shared_informative, len(shared_uninformative)))
            else:
                informative_subset = set()#if there are more shared uninformative than shared informative we will not subset
            # print(informative_subset)
            # exit(0)
        else:
            informative_subset = set()

        # print(informative_subset)
        self.informative_subset = informative_subset

        snp_gtypes_set = set(snp_gtypes)
        snp_gtypes_set = sorted(snp_gtypes_set)

        cellsnp_gtypes_set = set(cellsnp_gtypes)
        cellsnp_gtypes_set = sorted(cellsnp_gtypes_set)

        #for i in range(0, len(snp_gtypes)):
        for i in range(0, len(snp_gtypes_set)):
            discordant = False
            # snp_data = snp_gtypes[i].split('_')
            # cellsnp_data = cellsnp_gtypes[i].split('_')
            snp_data = snp_gtypes_set[i].split('_')
            cellsnp_data = cellsnp_gtypes_set[i].split('_')

            # the below will no longer work due to differing length of input strings
            # snp_alleles = [snp_gtypes[i][-3], snp_gtypes[i][-1]]
            # cellsnp_alleles = [cellsnp_gtypes[i][-3], cellsnp_gtypes[i][-1]]

            snp_alleles = [snp_data[4][0], snp_data[4][2]]
            cellsnp_alleles = [cellsnp_data[4][0], cellsnp_data[4][2]]

            snp_alleles_set = set(snp_alleles)
            cellsnp_alleles_set = set(cellsnp_alleles)
            
            snp_var = ('_').join(snp_data[0:4])
            cellsnp_var = ('_').join(cellsnp_data[0:4])

            if not cellsnp_var == snp_var:
                print("Error with strict discordance calculations: " + snp_gtypes[i] + " " + cellsnp_gtypes[i])
                exit(1)
            else:
                for allele in cellsnp_alleles_set:
                    if not allele in snp_alleles_set:#if a cellSNP allele is found that is not in the array data this is discordant
                        discordant = True
            
            if discordant == True:
                true_discordant+=1
                discordant_vars.append(cellsnp_var)
                if snp_var in self.uninformative_sites:
                    true_discordant_uninformative+=1
                    true_discordant_uninformative_ids.append(snp_var)
                elif snp_var in self.informative_sites:
                    true_discordant_informative+=1
                    true_discordant_informative_ids.append(snp_var)
            else:
                relaxed_concordant+=1
                concordant_vars.append(cellsnp_var)
                if snp_var in self.uninformative_sites:
                    relaxed_concordant_uninformative+=1
                    relaxed_concordant_uninformative_ids.append(snp_var)
                elif snp_var in self.informative_sites:
                    relaxed_concordant_informative+=1
                    relaxed_concordant_informative_ids.append(snp_var)
            if len(shared_uninformative) > 0:
                if snp_var in informative_subset:
                    if discordant == True:
                        subset_informative_discordant+=1
                    else:
                        subset_informative_concordant+=1
            #    true_discordant_count, relaxed_concordant_count, relaxed_concordant_informative_count, relaxed_concordant_uninformative_count, true_discordant_informative_count, true_discordant_uninformative_count, subset_informative_sites_concordant_count, subset_informative_sites_discordant_count 
        cell_vars2 = cell_vars.set_index('ids') 
        disc = pd.DataFrame(set(cell_vars2.loc[discordant_vars]['combo']),columns=['combo_x'])
        df_cd = pd.merge(cell_vars, expected_vars, how='inner', on = 'pos')
        disc2= pd.merge(disc, df_cd, how='inner', on = 'combo_x')
        disc2['expected_retrieved'] = disc2['0_x']+'::'+disc2['0_y']
        disc_sites_string = ';'.join(disc2['expected_retrieved'])
        return true_discordant, relaxed_concordant, relaxed_concordant_informative_ids, relaxed_concordant_uninformative_ids, true_discordant_informative_ids, true_discordant_uninformative_ids, discordant_vars, concordant_vars, disc_sites_string
    
    def read_concordance_calc(self,expected_vars,cell_vars):
            
            # This is a wrapper to add up the discordant reads in the cellsnp file.
            
            # expected genotype 0/0
            expected_hom_ref = expected_vars[expected_vars['vars'] == '0/0']
            hom_ref_sites = set(expected_hom_ref['ids'])
            cell_vars2 = cell_vars[cell_vars['ids'].isin(hom_ref_sites)]
            ad_hom_ref = cell_vars2['AD'].sum()
            oth_hom_ref = cell_vars2['OTH'].sum() 
            discordant_hom_ref = ad_hom_ref + oth_hom_ref
                    
            # expected genotype 0/1 or 1/0
            hets = ['0/1', '1/0']
            expected_het = expected_vars[expected_vars['vars'].isin(hets)]
            het_sites = set(expected_het['ids'])
            cell_vars3 = cell_vars[cell_vars['ids'].isin(het_sites)]
            discordant_het = cell_vars3['OTH'].sum()
            
            # expected genotype 1/1
            expected_hom_alt = expected_vars[expected_vars['vars'] == '1/1']
            hom_alt_sites = set(expected_hom_alt['ids'])
            cell_vars4 = cell_vars[cell_vars['ids'].isin(hom_alt_sites)]
            # DP + OTH - AD
            ad_hom_alt = cell_vars4['AD'].sum()
            dp_hom_alt = cell_vars4['DP'].sum()
            oth_hom_alt = cell_vars4['OTH'].sum()
            discordant_hom_alt = (dp_hom_alt + oth_hom_alt) - ad_hom_alt
            
            # Total analysis
            discordant_reads =  discordant_hom_ref + discordant_het + discordant_hom_alt
            total_dp = cell_vars['DP'].sum()
            total_oth = cell_vars['OTH'].sum()
            total_reads = total_dp + total_oth
            
            return total_reads,total_dp,total_oth,discordant_reads

    def read_condordance(self, expected_vars, cell_vars,discordant_vars, concordant_vars):
        '''
        get read level concordance using DP, AD and OTH format fields
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="total counts for ALT and REF">
        ##FORMAT=<ID=AD,Number=1,Type=Integer,Description="total counts for ALT">
        ##FORMAT=<ID=OTH,Number=1,Type=Integer,Description="total counts for other bases from REF and ALT">
        '''
        if not len(expected_vars) == len(cell_vars):
            print("length mismatch between expected vars and cell vars")
            exit(1)

        total_sites = len(expected_vars)
        #add cols for DP, AD< OTH
        cell_vars['DP'] = cell_vars[0].str.split("_").str[5].astype(int)
        cell_vars['AD'] = cell_vars[0].str.split("_").str[6].astype(int)
        cell_vars['OTH'] = cell_vars[0].str.split("_").str[7].astype(int)
        

        
        # Total
        total_reads,total_dp,total_oth,discordant_reads = self.read_concordance_calc(expected_vars,cell_vars)
        
        # uninformative
        cell_vars_uninformative = cell_vars[cell_vars['ids'].isin(self.uninformative_sites)]
        total_reads_uninformative,total_dp_uninformative,total_oth_uninformative,discordant_reads_uninformative = self.read_concordance_calc(expected_vars,cell_vars_uninformative)
        
        # informative 
        cell_vars_informative = cell_vars[cell_vars['ids'].isin(self.informative_sites)]
        total_reads_informative,total_dp_informative,total_oth_informative,discordant_reads_informative = self.read_concordance_calc(expected_vars,cell_vars_informative)
        
        # Split into concordant and discordant sites
        # concordant
        concordant_sites = cell_vars[cell_vars['ids'].isin(set(concordant_vars))]
        total_reads_for_concordant_sites,total_dp_for_concordant_sites,total_oth_for_concordant_sites,discordant_reads_for_concordant_sites = self.read_concordance_calc(expected_vars,concordant_sites)
        
        # discordant
        discordant_sites = cell_vars[cell_vars['ids'].isin(set(discordant_vars))]
        total_reads_for_discconcordant_sites,total_dp_for_discconcordant_sites,total_oth_for_discconcordant_sites,discordant_reads_for_discconcordant_sites = self.read_concordance_calc(expected_vars,discordant_sites)
              
        # Subset analysis
        cell_vars_informative_subset = cell_vars[cell_vars['ids'].isin(self.informative_subset)]
        total_reads_informative_subset,total_dp_informative_subset,total_oth_informative_subset,discordant_reads_informative_subset = self.read_concordance_calc(expected_vars,cell_vars_informative_subset)
  
        return total_sites, \
            self.informative_covered, \
                self.uninformative_covered, \
                    total_reads, \
                        discordant_reads, \
                            total_reads_informative, \
                                discordant_reads_informative, \
                                    total_reads_uninformative, \
                                        discordant_reads_uninformative, \
                                            total_reads_informative_subset, \
                                                    discordant_reads_informative_subset, \
                                                        total_reads_for_concordant_sites, \
                                                            discordant_reads_for_concordant_sites, \
                                                                total_reads_for_discconcordant_sites, \
                                                                    discordant_reads_for_discconcordant_sites
                                                                    
                                                    
    
    def get_discordance(self,expected_vars2,cell_vars2):
        Concordant_Sites = set(cell_vars2['combo']).intersection(set(expected_vars2['combo']))
        Discordant_sites = set(cell_vars2['combo'])-set(expected_vars2['combo'])  
        disc = pd.DataFrame(Discordant_sites,columns=['combo_x'])
        df_cd = pd.merge(cell_vars2, expected_vars2, how='inner', on = 'pos')
        disc2= pd.merge(disc, df_cd, how='inner', on = 'combo_x')
        disc2['expected_retrieved'] = disc2['0_x']+'::'+disc2['0_y']
        disc_sites = ';'.join(disc2['expected_retrieved'])
        return Concordant_Sites,Discordant_sites,disc_sites
    
    
    def retrieve_concordant_discordant_sites(self,expected_vars_norm,cell_vars):
        # This function has been inspired by Hails Concordance implementations, however hail has a pitfall that it performs a lot of other stuff under hood and requires intermediate sorting operations.
        # Since the single cell calculations requires concordance calculations per cell this becomes very computationally heavy on Hail, hence we have implemented concordance calculations here as part of the pipeline.
        # Author: M.Ozols
        
        cell_vars_norm = self.norm_genotypes(cell_vars)

        if len(cell_vars_norm) > 0:
            Total_Overlapping_sites = set(expected_vars_norm['ids']).intersection(set(cell_vars_norm['ids']))
            expected_vars2 = expected_vars_norm[expected_vars_norm['ids'].isin(Total_Overlapping_sites)]
            cell_vars2 = cell_vars_norm[cell_vars_norm['ids'].isin(Total_Overlapping_sites)]
            
            # Find exact discordant sites
            Concordant_Sites, Discordant_sites, _ = self.get_discordance(expected_vars2, cell_vars2)
            #find truly discordant sites
            true_discordant_count, relaxed_concordant_count, relaxed_concordant_informative_count, relaxed_concordant_uninformative_count, true_discordant_informative_count, true_discordant_uninformative_count, discordant_vars, concordant_vars, disc_sites_string = self.get_strict_discordance(expected_vars2, cell_vars2)
            #find discordant reads
            total_sites, \
                informative_sites, \
                    uninformative_sites, \
                        total_reads, \
                            discordant_reads, \
                                total_reads_informative, \
                                    discordant_reads_informative, \
                                        total_reads_uninformative, \
                                            discordant_reads_uninformative, \
                                                    total_reads_informative_subset, \
                                                        discordant_reads_informative_subset, \
                                                            total_reads_for_concordant_sites, \
                                                                discordant_reads_for_concordant_sites, \
                                                                    total_reads_for_discconcordant_sites, \
                                                                        discordant_reads_for_discconcordant_sites   = self.read_condordance(expected_vars2, cell_vars2, discordant_vars, concordant_vars)
                                                                        
            discordant_read_fraction_in_concordant_sites = f"{discordant_reads_for_concordant_sites}/{total_reads_for_concordant_sites}"
            discordant_read_fraction_in_discordant_sites = f"{discordant_reads_for_discconcordant_sites}/{total_reads_for_discconcordant_sites}"
            discordant_reads_uninformative_fraction = f"{discordant_reads_uninformative}/{total_reads_uninformative}"
            discordant_reads_informative_fraction = f"{discordant_reads_informative}/{total_reads_informative}"
            
            # sanity checks
            if total_reads != total_reads_for_concordant_sites+total_reads_for_discconcordant_sites:
                print("Error: total reads dont add up ")
                exit(1)
            if discordant_reads != discordant_reads_for_concordant_sites+discordant_reads_for_discconcordant_sites:
                print("Error: discordant reads dont add up ")
                exit(1)
            

        else:
            Total_Overlapping_sites = set()
            Concordant_Sites = set()
            Discordant_sites = set()
            disc_sites = ''
            true_discordant_count = 0
            relaxed_concordant_count = 0
            total_sites = 0
    
            discordant_reads = 0

        return Concordant_Sites, \
            Discordant_sites, \
            Total_Overlapping_sites, \
            disc_sites_string, \
            cell_vars_norm, \
            true_discordant_count, \
            relaxed_concordant_count, \
            relaxed_concordant_informative_count, \
            relaxed_concordant_uninformative_count, \
            true_discordant_informative_count, \
            true_discordant_uninformative_count, \
            total_sites, \
            informative_sites, \
            uninformative_sites, \
            total_reads, \
            total_reads_informative, \
            total_reads_uninformative, \
            discordant_reads, \
            discordant_reads_informative, \
            discordant_reads_uninformative,  \
            discordant_vars, \
            concordant_vars, \
            discordant_read_fraction_in_concordant_sites, \
            discordant_read_fraction_in_discordant_sites, discordant_reads_uninformative_fraction, discordant_reads_informative_fraction
    

    def set_results(self,to_set,id):
        # Recod to disk to save the loading mmeory time.
        with open(f'tmp_{id}.pkl', 'wb') as f:
            pickle.dump(to_set, f)
        self.record_dict[id]=f'tmp_{id}.pkl'
        return 
    
    # def append_results_cell_concordances(self,result):
    def append_results_cell_concordances(self,result,cell_concordance_table,other_donor_concordances,other_donor_concordance_table):
        other_donor_concordance_table = other_donor_concordance_table + other_donor_concordances
        count=result['count']
        try:
            percent_concordant = result['Nr_Concordant']/(result['Nr_Discordant']+result['Nr_Concordant'])*100
        except:
            percent_concordant = 0
        
        try:
            percent_discordant = result['Nr_Discordant']/(result['Nr_Discordant']+result['Nr_Concordant'])*100
        except:
            percent_discordant = 0

        try:
            percent_relaxed_concordant = result['Nr_Relaxed_concordant']/(result['Nr_Relaxed_concordant']+result['true_discordant_count'])*100
        except:
            percent_relaxed_concordant = 0
        
        try:
            percent_strict_discordant = result['true_discordant_count']/(result['Nr_Relaxed_concordant']+result['true_discordant_count'])*100
        except:
            percent_strict_discordant = 0

        try:
            read_discordance = result['discordant_reads']/result['total_sites']
        except:
            read_discordance = 0

        cohort = 'UNKNOWN'
        donor_split = result['donor_gt_match'].split("_")
        if (len(donor_split) == 2) and (donor_split[0] == donor_split[1]):
            cohort = 'UKB'
        elif (len(donor_split) == 3) and (len(donor_split[0]) == 14):
            cohort = 'ELGH'

        same_as_asigned_donor = result['donor_gt_match'] in result['Donor_With_Highest_Concordance']
        if not same_as_asigned_donor:
            same_as_asigned_donor = result['donor_gt_match'] in result['Donor_With_Lowest_DisConcordance']
            
        cell_concordance_table[f"{result['cell1']} --- {result['donor_gt_match']}"] = {'GT 1':result['cell1'],
                                                                'GT 2':result['donor_gt_match'],
                                                                'cohort': cohort,
                                                                
                                                                'Nr_Concordant':result['Nr_Concordant'],
                                                                'Nr_Discordant':result['Nr_Discordant'],
                                                                'Nr_Relaxed_concordant':result['Nr_Relaxed_concordant'],
                                                                'Nr_strict_discordant':result['true_discordant_count'],
                                                                'Percent Concordant':percent_concordant,
                                                                'Percent Discordant':percent_discordant,
                                                                'Percent_relaxed_concordant': percent_relaxed_concordant,
                                                                'Percent_strict_discordant': percent_strict_discordant,
                                                                'Nr_concordant_informative': len(result['relaxed_concordant_informative_count']),
                                                                'Nr_concordant_uninformative': len(result['relaxed_concordant_uninformative_count']),
                                                                'Nr_discordant_informative': len(result['true_discordant_informative_count']),
                                                                'Nr_discordant_uninformative': len(result['true_discordant_uninformative_count']),
                                                                'NrTotal_Overlapping_sites_between_two_genotypes':result['Nr_Total_Overlapping_sites'],
                                                                'Nr_donor_distinct_sites_within_pool_individuals':result['Nr_donor_distinct_sites'],
                                                                'Number_of_sites_that_are_donor_concordant_and_exclusive':result['Number_of_sites_that_are_donor_concordant_and_exclusive'],
                                                                'Total_sites': result['total_sites'],
                                                                'Total_informative_sites': result['informative_sites'],
                                                                'Total_uninformative_sites': result['uninformative_sites'],
                                                                'Total_reads': result['total_reads'],
                                                                'Total_reads_informative': result['total_reads_informative'],
                                                                'Total_reads_uninformative': result['total_reads_uninformative'],
                                                                'Discordant_reads': result['discordant_reads'],
                                                                'Discordant_reads_informtive': result['discordant_reads_informative'],
                                                                'Discordant_reads_uninformtive': result['discordant_reads_uninformative'],
                                                                'Discordant_reads_by_n_sites': read_discordance,
                                                                
                                                                'Discordant_sites_in_pool': len(result['Discordant_sites_in_pool']),
                                                                'Lowest_Disconcordance_value_in_all_donors':result['Lowest_Disconcordance_value_in_all_donors'],
                                                                'Donor_With_Lowest_DisConcordance':result['Donor_With_Lowest_DisConcordance'],
                                                                'Concordant_Site_Identities':result['Concordant_Site_Identities'],
                                                                'Donor_With_Highest_Concordance':result['Donor_With_Highest_Concordance'],
                                                                'Highest_Concordance_value_in_all_donors':result['Highest_Concordance_value_in_all_donors'],
                                                                'same_as_asigned_donor':same_as_asigned_donor,
                                                                'Total_sites_other_donor (if same_as_asigned_donor=False)':result['Total_sites_other_donor'],
                                                                'Total_reads_other_donor (if same_as_asigned_donor=False)':result['Total_reads_other_donor'],
                                                                'total_discordant_sites_that_are_concordant_with_other_donors_in_pool':result['total_discordant_sites_that_are_concordant_with_other_donors_in_pool'],
                                                                'discordant_read_fraction_in_concordant_site':result['discordant_read_fraction_in_concordant_sites'], 
                                                                'discordant_read_fraction_in_discordant_sites':result['discordant_read_fraction_in_discordant_sites'], 
                                                                'Discordant_Site_Identities':result['discordant_sites'],
                                                                } 
        
        return [cell_concordance_table,other_donor_concordance_table]
    
    # def combine_written_files(self):#this one is for concordance class
    #     to_export = self.cell_concordance_table
    #     for val1 in self.record_dict.values():
    #         # here remove the int files.
    #         print(f"merging temp file: {val1}")
    #         with open(val1, 'rb') as f:
    #             loaded_dict = pickle.load(f)
    #             for k1 in loaded_dict.keys():
    #                 to_export[k1]=loaded_dict[k1]
    #         os.remove(val1)
    #     return to_export
    
    
    def combine_written_lists(self,exclusive_donor_variants,record_dict):#this is for VCF loader class
        to_export = exclusive_donor_variants
        for val1 in record_dict.values():
            # here remove the int files.
            print(f"merging temp file: {val1}")
            with open(val1, 'rb') as f:
                loaded_dict = pickle.load(f)
                self.other_donor_comp = self.other_donor_comp+ loaded_dict
            os.remove(val1)
        return self.other_donor_comp    
    
    def combine_written_files(self,exclusive_donor_variants,record_dict):#this is for VCF loader class
        to_export = exclusive_donor_variants
        for val1 in record_dict.values():
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
   
    def set_results(self,to_set,id):
        # Recod to disk to save the loading mmeory time.
        with open(f'tmp_{id}.pkl', 'wb') as f:
            pickle.dump(to_set, f)
        self.record_dict[id]=f'tmp_{id}.pkl' 
    
    def analyse_donor(self,Cells_to_keep_pre,donor_gt_match,donor_gt_match_cohort,vars_per_donor_gt,donor_cohorts,count,all_donor_data,expected_vars_norm):
        donor_concordance_table = {}
        other_donor_concordance_table = []
        for cell1 in Cells_to_keep_pre:
            count+=1

            cell_vars = exclusive_cell_variants[cell1]
            result1, other_donor_concordances = self.concordance_table_production(expected_vars_norm,cell_vars,cell1,donor_gt_match,donor_gt_match_cohort, vars_per_donor_gt, donor_cohorts, count,all_donor_data)
            cell_concordance_table,other_donor_concordance_table = self.append_results_cell_concordances(result1,donor_concordance_table,other_donor_concordances,other_donor_concordance_table)
            if count>300:
                break
            # here we should write these independently to the files
            if (count % 50 == 0):
                self.set_results(other_donor_concordance_table,f"{count}--{donor_gt_match}")
                other_donor_concordance_table = []
                

        self.set_results(other_donor_concordance_table,f"{count}--{donor_gt_match}")
        output2 = self.combine_written_lists(self.other_donor_comp,self.record_dict)
        pd.DataFrame(output2).sort_values(by=['cell']).to_csv(f'{donor_gt_match}--each_cells_comparison_with_other_donor.tsv',sep='\t',index=False)
        del output2
        return donor_concordance_table

    def combine_concordances(self,result):
        
        self.cell_concordance_table = {**self.cell_concordance_table, **result}
        

    def conc_table(self):
        donor_assignments_table=self.donor_assignments_table
        cell_assignments_table=self.cell_assignments_table
        exclusive_don_variants=self.exclusive_don_variants
        exclusive_cell_variants= self.exclusive_cell_variants
        donor_list = exclusive_don_variants.keys()
        pool = mp.Pool(cpus)
        count = 0
        
        
        #create a list of variants that are on each donor genotype file
        vars_per_donor_gt = {}
        for don_id in donor_list:
            donor_gt_vars = list(exclusive_don_variants[don_id])
            donor_gt_vars = pd.DataFrame(donor_gt_vars)
            donor_gt_vars = self.norm_genotypes(donor_gt_vars)
            donor_gt_vars = donor_gt_vars[donor_gt_vars['vars'] != '0/0']
            donor_gt_varids = list(donor_gt_vars['ids'])
            vars_per_donor_gt[don_id] = donor_gt_varids
        
        #work out what cohort each donor belongs to
        donor_cohorts = {}
        for don_id in donor_list:
            cohort = 'UNKNOWN'
            donor_split = don_id.split("_")
            if (len(donor_split) == 2) and (donor_split[0] == donor_split[1]):
                cohort = 'UKB'
            elif (len(donor_split) == 3) and (len(donor_split[0]) == 14):
                cohort = 'ELGH'
            donor_cohorts[don_id] = cohort

        all_donor_data={}
        # here we calvculate all the expected donor datasets
        for row1 in exclusive_don_variants.keys():
            # donor_in_question = row1['donor_query']
            donor_gt_match = row1
            expected_vars_of_other_donor = self.exclusive_don_variants[donor_gt_match]
            expected_vars_norm_of_other_donor = self.norm_genotypes(expected_vars_of_other_donor)
            all_donor_data[donor_gt_match]=expected_vars_norm_of_other_donor

        for i,row1 in donor_assignments_table.iterrows():
            donor_in_question = row1['donor_query']
            donor_gt_match = row1['donor_gt']
            # if i>4:
            #     continue
            if (donor_gt_match=='NONE'):
                continue
            try:
                donor_gt_match_cohort = donor_cohorts[donor_gt_match]
            except:
                continue
            Cells_to_keep_pre = list(set(cell_assignments_table.loc[cell_assignments_table['donor_id']==donor_in_question,'cell']))
            expected_vars = exclusive_don_variants[donor_gt_match]
            expected_vars_norm = self.norm_genotypes(expected_vars)
            try:
                # Now we subset this down to each of the uniqie variants per donor and check which of the concordant sites are exclusive to donor.
                dds = self.donor_distinct_sites[donor_gt_match]
            except:
                continue
            if cpus==1:
                result = self.analyse_donor(Cells_to_keep_pre,donor_gt_match,donor_gt_match_cohort,vars_per_donor_gt,donor_cohorts,count,all_donor_data,expected_vars_norm)
                self.combine_concordances(result)
            else:
                pool.apply_async(self.analyse_donor, args=([Cells_to_keep_pre,donor_gt_match,donor_gt_match_cohort,vars_per_donor_gt,donor_cohorts,count,all_donor_data,expected_vars_norm]),callback=self.combine_concordances)          
                
        pool.close()
        pool.join()
        
        # output = self.combine_written_files(self.cell_concordance_table,self.record_dict)
        
        return self.cell_concordance_table

    
    def concordance_table_production(self,expected_vars_norm,cell_vars,cell1,donor_gt_match, donor_gt_match_cohort, vars_per_donor_gt, donor_cohorts, count,all_donor_data):

        Concordant_Sites, \
            Discordant_sites, \
            Total_Overlapping_sites, \
            discordant_sites, \
            cell_vars_norm, \
            true_discordant_count, \
            relaxed_concordant_count, \
            relaxed_concordant_informative_count, \
            relaxed_concordant_uninformative_count, \
            true_discordant_informative_count, \
            true_discordant_uninformative_count, \
            total_sites, \
            informative_sites, \
            uninformative_sites, \
            total_reads, \
            total_reads_informative, \
            total_reads_uninformative, \
            discordant_reads, \
            discordant_reads_informative, \
            discordant_reads_uninformative, \
            discordant_vars, \
            concordant_vars, \
            discordant_read_fraction_in_concordant_sites, \
            discordant_read_fraction_in_discordant_sites, \
            discordant_reads_uninformative_fraction, \
            discordant_reads_informative_fraction  = self.retrieve_concordant_discordant_sites(expected_vars_norm,cell_vars)
            
        total_concordant_sites = len(Concordant_Sites) + relaxed_concordant_count
        dds = self.donor_distinct_sites[donor_gt_match]
        Nr_donor_distinct_sites = len(dds)
        Nr_Concordant = len(Concordant_Sites)
        Nr_Relaxed_concordant = Nr_Concordant + relaxed_concordant_count
        Nr_Discordant = len(Discordant_sites)
        Nr_Total_Overlapping_sites = len(Total_Overlapping_sites)
        Number_of_sites_that_are_donor_concordant_and_exclusive = len(set(dds).intersection(set(Discordant_sites)))
        Number_of_sites_in_cellsnp_but_not_in_reference = set(cell_vars_norm['pos'])-set(expected_vars_norm['pos'])
        #Quantify donor variation in other donors
        discordant_vars_in_pool = []
        donor_table_of_concordances = []
        total_discordant_sites_that_are_concordant_with_other_donors_in_pool = set()
        informative__total_discordant_sites_that_are_concordant_with_other_donors_in_pool = set()
        total_cordant_sites_that_are_concordant_with_other_donors_in_pool = set()
        for donor in vars_per_donor_gt:
          
            expected_vars_norm_of_other_donor = all_donor_data[donor]
            
            Concordant_Sites_otherDonor, \
                Discordant_sites_otherDonor, \
                Total_Overlapping_sites_otherDonor, \
                discordant_sites_otherDonor, \
                cell_vars_norm_otherDonor, \
                true_discordant_count_otherDonor, \
                relaxed_concordant_count_otherDonor, \
                relaxed_concordant_informative_count_otherDonor, \
                relaxed_concordant_uninformative_count_otherDonor, \
                true_discordant_informative_count_otherDonor, \
                true_discordant_uninformative_count_otherDonor, \
                total_sites_otherDonor, \
                informative_sites_otherDonor, \
                uninformative_sites_otherDonor, \
                total_reads_otherDonor, \
                total_reads_informative_otherDonor, \
                total_reads_uninformative_otherDonor, \
                discordant_reads_otherDonor, \
                discordant_reads_informative_otherDonor, \
                discordant_reads_uninformative_otherDonor, \
                discordant_vars_otherDonor, \
                concordant_vars_otherDonor, \
                discordant_read_fraction_in_concordant_sites_otherDonor, \
                discordant_read_fraction_in_discordant_sites_otherDonor, \
                discordant_reads_uninformative_fraction_otherDonor, \
                discordant_reads_informative_fraction_otherDonor  = self.retrieve_concordant_discordant_sites(expected_vars_norm_of_other_donor,cell_vars)
            
            # here we also need to know :
            #   how many reads of the desired donor discordant sites could be yielded
   

            total_concordant_sites_otherDonor = relaxed_concordant_count_otherDonor
            concordant_percent_in_other_donor= total_concordant_sites_otherDonor/total_sites_otherDonor*100
            discordant_percent_in_other_donor= true_discordant_count_otherDonor/total_sites_otherDonor*100
            DonorDiscordant_Sites_that_are_atributed_to_other_donor = set(discordant_vars).intersection(set(concordant_vars_otherDonor))
            Informative__DonorDiscordant_Sites_that_are_atributed_to_other_donor = set(true_discordant_informative_count).intersection(set(relaxed_concordant_informative_count_otherDonor))
            DonorCordant_Sites_that_are_atributed_to_other_donor = set(concordant_vars).intersection(set(concordant_vars_otherDonor))

            # We now count the concordant reads  that may contribute to  particular cell at this cell.
            # to do this we take the discordant sites that have been deamed to be concordant with the other donor and quantify the reads thta are concordant.
            Total_Overlapping_sites = set(DonorDiscordant_Sites_that_are_atributed_to_other_donor)
            expected_vars2 = expected_vars_norm_of_other_donor[expected_vars_norm_of_other_donor['ids'].isin(Total_Overlapping_sites)]
            cell_vars2 = cell_vars_norm[cell_vars_norm['ids'].isin(Total_Overlapping_sites)]        
            cell_vars2['DP'] = cell_vars2[0].str.split("_").str[5].astype(int)
            cell_vars2['AD'] = cell_vars2[0].str.split("_").str[6].astype(int)
            cell_vars2['OTH'] = cell_vars2[0].str.split("_").str[7].astype(int)
            
            total_reads_for_discordant_sites_that_are_concordant_with_other_donor,total_dp_for_discordant_sites_that_are_concordant_with_other_donor,total_oth_for_discordant_sites_that_are_concordant_with_other_donor,discordant_reads_for_discordant_sites_that_are_concordant_with_other_donor = self.read_concordance_calc(expected_vars2,cell_vars2)
            concordant_reads_for_discordant_sites_that_are_concordant_with_other_donor = total_reads_for_discordant_sites_that_are_concordant_with_other_donor - discordant_reads_for_discordant_sites_that_are_concordant_with_other_donor
            
            try:
                donor_cohort = donor_cohorts[donor]
                donor_vars = vars_per_donor_gt[donor]
            except:
                continue            
            if not donor == donor_gt_match:
                # We want to kow how many of these discordant site

                total_discordant_sites_that_are_concordant_with_other_donors_in_pool = total_discordant_sites_that_are_concordant_with_other_donors_in_pool.union(set(DonorDiscordant_Sites_that_are_atributed_to_other_donor))
                # to get the total reads that can be atributed to the other donor i have to check if site is already covered in the total_discordant_sites_that_are_concordant_with_other_donors_in_pool.
                # the ones that havent, i have to add the reads up for them.
                informative__total_discordant_sites_that_are_concordant_with_other_donors_in_pool = informative__total_discordant_sites_that_are_concordant_with_other_donors_in_pool.union(set(Informative__DonorDiscordant_Sites_that_are_atributed_to_other_donor))
                
                total_cordant_sites_that_are_concordant_with_other_donors_in_pool = total_cordant_sites_that_are_concordant_with_other_donors_in_pool.union(set(DonorCordant_Sites_that_are_atributed_to_other_donor))
                

                common_vars = list(set(discordant_vars) & set(donor_vars))
                common_var_count = str(len(common_vars))
                donor_cohort_common = donor + ":" + donor_cohort + ":" + common_var_count
                discordant_vars_in_pool.append(donor_cohort_common)
                
                # Here we want to calculate the number of discordant sites in other donors and see if in terms of concordance the same donor is picked as per GT assignment.
                # We do this to investigate the potential of a cell coming from this other donor.
            
            donor_table_of_concordances.append({'donor':donor, 'cell':cell1, 'donor_cohort':donor_cohort, \
                                                'gt matched donor':donor == donor_gt_match, \
                                                'DonorCordant_Sites_that_are_atributed_to_other_donor':len(DonorCordant_Sites_that_are_atributed_to_other_donor), \
                                                'DonorCordant_Sites_that_are_atributed_to_other_donor/total':f"{len(DonorCordant_Sites_that_are_atributed_to_other_donor)}/{len(concordant_vars)}", \
                                                'DonorDiscordant_Sites_that_are_atributed_to_other_donor':len(DonorDiscordant_Sites_that_are_atributed_to_other_donor), \
                                                'DonorDiscordant_Sites_that_are_atributed_to_other_donor/total':f"{len(DonorDiscordant_Sites_that_are_atributed_to_other_donor)}/{len(discordant_vars)}", \
                                                'concordant_percent_in_other_donor':concordant_percent_in_other_donor, \
                                                'discordant_percent_in_other_donor':discordant_percent_in_other_donor, \
                                                'discordant_reads_otherDonor':discordant_reads_otherDonor, \
                                                'discordant_sites_otherDonor':len(discordant_vars_otherDonor), \
                                                'concordant_sites_otherDonor':len(concordant_vars_otherDonor), \
                                                'total_sites_otherDonor':total_sites_otherDonor, \
                                                'discordant_reads_otherDonor':discordant_reads_otherDonor, \
                                                'total_reads_otherDonor':total_reads_otherDonor, \
                                                'discordant_read_fraction_in_concordant_sites_otherDonor':discordant_read_fraction_in_concordant_sites_otherDonor, \
                                                'discordant_read_fraction_in_discordant_sites_otherDonor':discordant_read_fraction_in_discordant_sites_otherDonor, \
                                                'concordant_reads_for_discordant_sites_that_are_concordant_with_other_donor':concordant_reads_for_discordant_sites_that_are_concordant_with_other_donor
                                                })
                
        discordant_vars_in_pool_str = (";").join(discordant_vars_in_pool)
        concordant_vars_in_pool_str = (";").join(concordant_vars)
        DF = pd.DataFrame(donor_table_of_concordances)
        
        Donor_With_Lowest_DisConcordance = ';'.join(DF[DF['discordant_percent_in_other_donor']==min(DF['discordant_percent_in_other_donor'])]['donor'].values)
        Lowest_Disconcordance_value_in_all_donors= DF[DF['discordant_percent_in_other_donor']==min(DF['discordant_percent_in_other_donor'])]['discordant_percent_in_other_donor'].values[0]
        
        Donor_With_Highest_Concordance = ';'.join(DF[DF['concordant_percent_in_other_donor']==max(DF['concordant_percent_in_other_donor'])]['donor'].values)
        Highest_Concordance_value_in_all_donors= DF[DF['concordant_percent_in_other_donor']==max(DF['concordant_percent_in_other_donor'])]['concordant_percent_in_other_donor'].values[0]
        Total_sites_other_donor = ';'.join(DF[DF['concordant_percent_in_other_donor']==max(DF['concordant_percent_in_other_donor'])]['total_sites_otherDonor'].astype(str).values)
        Total_reads_other_donor = ';'.join(DF[DF['concordant_percent_in_other_donor']==max(DF['concordant_percent_in_other_donor'])]['total_reads_otherDonor'].astype(str).values)
                                                                     
        return [{
            'cell1':cell1,
            'donor_gt_match':donor_gt_match,
            'Nr_Concordant':Nr_Concordant,
            'Nr_Discordant':Nr_Discordant,
            'Nr_Relaxed_concordant':Nr_Relaxed_concordant,
            'true_discordant_count':true_discordant_count,
            'relaxed_concordant_informative_count':relaxed_concordant_informative_count,
            'relaxed_concordant_uninformative_count':relaxed_concordant_uninformative_count,
            'true_discordant_informative_count':true_discordant_informative_count,
            'true_discordant_uninformative_count':true_discordant_uninformative_count,
            'Nr_Total_Overlapping_sites':Nr_Total_Overlapping_sites,
            'Number_of_sites_that_are_donor_concordant_and_exclusive':Number_of_sites_that_are_donor_concordant_and_exclusive,
            'Nr_donor_distinct_sites':Nr_donor_distinct_sites,
            'count':count,
            'discordant_sites':discordant_sites,
            'total_sites':total_sites,
            'informative_sites':informative_sites,
            'uninformative_sites':uninformative_sites,
            'total_reads':total_reads,
            'total_reads_informative':total_reads_informative,
            'total_reads_uninformative':total_reads_uninformative,
            'discordant_reads':discordant_reads,
            'discordant_reads_informative':discordant_reads_informative,
            'discordant_reads_uninformative':discordant_reads_uninformative,
            'Discordant_sites_in_pool': discordant_vars,
            'Lowest_Disconcordance_value_in_all_donors':Lowest_Disconcordance_value_in_all_donors,
            'Donor_With_Lowest_DisConcordance':Donor_With_Lowest_DisConcordance,
            'Concordant_Site_Identities':concordant_vars_in_pool_str,
            'Donor_With_Highest_Concordance':Donor_With_Highest_Concordance,
            'Highest_Concordance_value_in_all_donors':Highest_Concordance_value_in_all_donors,
            'Total_sites_other_donor':Total_sites_other_donor,
            'Total_reads_other_donor':Total_reads_other_donor,
            'total_discordant_sites_that_are_concordant_with_other_donors_in_pool':f"{len(total_discordant_sites_that_are_concordant_with_other_donors_in_pool)}/{len(discordant_vars)}",
            'informative__total_discordant_sites_that_are_concordant_with_other_donors_in_pool':f"{len(total_discordant_sites_that_are_concordant_with_other_donors_in_pool)}/{len(true_discordant_informative_count)}",
            'discordant_read_fraction_in_concordant_sites':discordant_read_fraction_in_concordant_sites, \
            'discordant_read_fraction_in_discordant_sites':discordant_read_fraction_in_discordant_sites
        }, donor_table_of_concordances]
    
    
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
        '''
        takes VCF lines and extracts all format fields for those where GT !='.'
        '''
        list_val = line.rstrip().split("\t") #[:5] #:8
        idx = find(list_val[8].split(':'),'GT')[0]#find index of GT field as GT will tell us what variants are called
        if remove_AG:
            if list_val[3] == 'A' and list_val[4] == 'G':#remove A>G
                pass
            elif list_val[3] == 'T' and list_val[4] == 'C':#also remove T>C
                pass
        
        if len(list_val[3]) > 1 or len(list_val[4]) > 1:
            # CURRENTLY DEALS ONLY WITH BIALELIC
            print(f'{idx} var not bialelic')
            pass
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

        return [obs_with_gt,list_val_with_gt,idx,list_val,count,format_list]#add format_list to the return value as we need this for the next step


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
        format_list = result[5]#list of required format fields
        #get indexes of required format fields (apart from GT which has already been taken care of)
        additional_field_idxs = []
        for fmt in format_list:
            if not fmt == 'GT':
                idx_addn = find(list_val[8].split(':'), fmt)[0]
                additional_field_idxs.append(idx_addn)
        # print(additional_field_idxs)
        # exit(0)

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
            #append any additional format fields to alleles
            if len(additional_field_idxs) > 0:
                for idx_addnl in additional_field_idxs:
                    fmt_val = list_val_with_gt[donor_loc_in_list].split(':')[idx_addnl]
                    alleles = alleles + '_' + fmt_val

            if not alleles.startswith('.'):
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

    def combine_written_files(self,exclusive_donor_variants,record_dict):#this is for VCF loader class
        to_export = exclusive_donor_variants
        for val1 in record_dict.values():
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
        
        output = self.combine_written_files(self.exclusive_donor_variants,self.record_dict)
        
        return output
    
"""Run CLI."""

def get_options():
    '''
    Get options from the command line
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--version', action='version', version='%(prog)s {version}'.format(version=__version__))
    parser.add_argument('--cpus', action='store', required=True, type=int)
    parser.add_argument('--cell_vcf', action='store', required=True)
    parser.add_argument('--cell_assignments', action='store', required=True)
    parser.add_argument('--donor_assignments', action='store', required=True)
    parser.add_argument('--gt_match_vcf', action='store', required=True)
    parser.add_argument('--expected_vcf', action='store', required=True)
    parser.add_argument('--informative_sites', action='store', required=True)
    parser.add_argument('--uninformative_sites', action='store', required=True)
    parser.add_argument('--outfile', action='store', required=True)
    parser.add_argument('--debug', action='store_true')
    parser.add_argument('--remove_AG', action='store_true')
    args = parser.parse_args()

    return args


def get_sites_from_tsv(sites_file):
    """
    get sites frm a tsv file where cols are chrom, pos, id, ref, alt
    assumes no multiallelics
    """
    sites = set()
    with open(sites_file, 'r') as f:
        lines = f.readlines()
        for l in lines:
            linedata = l.split('\t')
            var = ('_').join([linedata[0], linedata[1], linedata[3], linedata[4]])
            sites.add(var)
    return sites
    

def find(lst, a):
    return [i for i, x in enumerate(lst) if x==a ]
def norm_genotypes(expected_vars):
    expected_vars = pd.DataFrame(expected_vars)
    split_str=expected_vars[0].str.split("_")
    expected_vars['ids'] = split_str.str[0]+'_'+split_str.str[1]+'_'+split_str.str[2]+'_'+split_str.str[3]
    expected_vars['pos'] = split_str.str[0]+'_'+split_str.str[1]
    expected_vars['vars'] = split_str.str[4]
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




if __name__ == "__main__":

    options = get_options()
    cpus = options.cpus
    outfile = options.outfile
    cell_vcf=options.cell_vcf
    donor_assignments=options.donor_assignments
    gt_match_vcf=options.gt_match_vcf
    expected_vcf=options.expected_vcf
    cell_assignments=options.cell_assignments
    informative_sites_file = options.informative_sites
    uninformative_sites_file = options.uninformative_sites

    informative_sites = get_sites_from_tsv(informative_sites_file)
    uninformative_sites = get_sites_from_tsv(uninformative_sites_file)

    exclusive_donor_variants = {} #This is where results are populated when mp process i used.
    curently_pushing =[] #this is a lock value to check if rhe curent field is updated so to avaid the race for update
    All_Results={}
    cell_concordance_table = {}

    donor_assignments_table = pd.read_csv(donor_assignments)
    cell_assignments_table = pd.read_csv(cell_assignments,sep='\t')
    remove_AG = options.remove_AG
    
    if options.debug:
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
    else:  
        print('---Loading genotype VCF----')   
        if (os.path.exists(gt_match_vcf)):
            loader2 = VCF_Loader(gt_match_vcf, biallelic_only=True,
                            sparse=False, format_list=['GT'])
            GT_Matched_variants = loader2.load_VCF_batch_paralel()
            del loader2
        else:
            GT_Matched_variants = {}
        
        with open(f'tmp_GT_Matched_variants.pkl', 'wb') as f:
            pickle.dump(GT_Matched_variants, f)
        
        print('---Loading cell VCF----')
        loader1 = VCF_Loader(cell_vcf, biallelic_only=True,
                            sparse=False, format_list=['GT', 'DP', 'AD', 'OTH'])
        exclusive_cell_variants = loader1.load_VCF_batch_paralel()
        del loader1
        with open(f'tmp_exclusive_cell_variants.pkl', 'wb') as f:
            pickle.dump(exclusive_cell_variants, f)

        print('---Loading expected VCF----')
        loader3 = VCF_Loader(expected_vcf, biallelic_only=True,
                        sparse=False, format_list=['GT'])
        GT_Expected_variants = loader3.load_VCF_batch_paralel()
        del loader3

        with open(f'tmp_GT_Expected_variants.pkl', 'wb') as f:
            pickle.dump(GT_Expected_variants, f)
 
        print('---Variant files loaded----')       
        
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

    cell_concordance_table = Concordances(donor_assignments_table,cell_assignments_table,exclusive_don_variants,exclusive_cell_variants,donor_distinct_sites, informative_sites, uninformative_sites).conc_table()
    

    result = pd.DataFrame(cell_concordance_table).T
    try:
        site_identities = result[['Concordant_Site_Identities','Discordant_Site_Identities']]
        result.drop(columns=['Concordant_Site_Identities','Discordant_Site_Identities'],inplace=True)
        site_identities.to_csv(f"site_identities_{outfile}",sep='\t')
    except:
        _='sample_hasnt_matched_any_gt --- most likely too little cells assigned'
    result.to_csv(outfile,sep='\t')
    
    print('Processing Done')