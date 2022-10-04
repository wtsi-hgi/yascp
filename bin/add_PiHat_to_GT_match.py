#!/usr/bin/env python


__date__ = '2022-10-20'
__author__ = 'M.Ozols'
__version__ = '0.0.1'

import argparse
import pandas as pd
import re
def main():
    """Run CLI."""
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
       '-ph', '--pihat_match',
        action='store',
        dest='pihat_match',
        required=True,
        help='PiHAT values for match'
    )

    parser.add_argument(
       '-c', '--cond',
        action='store',
        dest='condition',
        required=True,
        help='colname for condition'
    )
    
    parser.add_argument(
       '-e', '--expected',
        action='store',
        dest='expected',
        required=True,
        help='expected sample IDs'
    )
    
    parser.add_argument(
       '-m', '--mapping_file',
        action='store',
        dest='mapping_file',
        required=False,
        default=None,
        help='Mapping between genotype/phenotype'
    )
    
    parser.add_argument(
       '-id', '--pool_id',
        action='store',
        dest='pool_id',
        required=False,
        default=None,
        help='Mapping between genotype/phenotype'
    )

    parser.add_argument(
       '-mt', '--match_table',
        action='store',
        dest='match_table',
        required=True,
        help='match_table'
    )

    options = parser.parse_args()
    Genome_PiHAT = pd.read_csv(options.pihat_match,sep='\s+')
    GT_Match_Table = pd.read_csv(options.match_table,sep=',')
    # Name = options.match_table.replace('stats_','').replace('_gt_donor_assignments.csv','')
    Name = options.pool_id
    Condition_Column = options.condition
    expected_ids = options.expected.split(',')
    mapping_file = options.mapping_file
    if mapping_file:
        bridge = pd.read_csv(mapping_file,sep='\t')
        all_maped_samples = []
        for s1 in expected_ids:
            # there are no dublicates atm for the S2 ids, but there are replicates for the genotype ids, this however not conceir us.
            mapping = bridge[bridge['s00046_id']==s1]
            if len(mapping)==0:
                s2 = re.sub('^0*','',s1)
                mapping = bridge[bridge['s00046_id']==s2]
            if len(mapping)==0:
                s2 = s1.split('_')[0]
                mapping = bridge[bridge['s00046_id']==s2]
                
            try:
                replacement = mapping['oragene_id'].values[0]
            except:
                replacement = s1
            all_maped_samples.append({'original':s1,'replacement':replacement})
        all_maped_samples2 = pd.DataFrame(all_maped_samples)
        # all_maped_samples2=all_maped_samples2.set_index('original')
    else:
        all_maped_samples2 = expected_ids
        
    if Condition_Column=='Expected':
        # Here we produce a new file that indicates all the best matches for the expected samples and their associated PiHat values.
        all_expected_id_best_match_PiHat = []
        for i,s1 in all_maped_samples2.iterrows():
            set1 =Genome_PiHAT[Genome_PiHAT['IID1'].str.contains(s1.replacement)]
            set2 =Genome_PiHAT[Genome_PiHAT['IID2'].str.contains(s1.replacement)]
            combo=pd.concat([set1,set2])
            try:
                max_pihat = combo['PI_HAT'].values.max()
                matching_donor_row = combo[combo['PI_HAT'] == max_pihat]
                if (matching_donor_row['IID1'].values[0]==s1.replacement):
                    donor_IID = matching_donor_row['IID2'].values[0]
                else:
                    donor_IID = matching_donor_row['IID1'].values[0]
                
            except:
                if (s1.original == s1.replacement):
                    donor_IID = 'likely that mapping file did not contain bridging info for this sample'
                else:
                    donor_IID = 'likely that sample mapping indicated is not in VCF file'
                max_pihat = None
            try:
                donor_nr =int(donor_IID.replace('donor',''))
            except:
                donor_nr = -1
            all_expected_id_best_match_PiHat.append({'expected_donor_id':s1.original,'Best_Matched_donor_IID':donor_IID,'PiHAT value':max_pihat, 'mapping_id':s1.replacement,'donor_nr':donor_nr})
                
        all_expected_data_PiHAT = pd.DataFrame(all_expected_id_best_match_PiHat)
        
        all_expected_data_PiHAT= all_expected_data_PiHAT.sort_values(by=['donor_nr','PiHAT value'])
        all_expected_data_PiHAT=all_expected_data_PiHAT.drop(columns=['donor_nr'])
        all_expected_data_PiHAT.to_csv(f'Max_PiHAT_For_Expected_{Name}.tsv',sep='\t',index=False)
        
    GT_Match_Table[f"PiHat: {Condition_Column}"]=None
    for i,row1 in GT_Match_Table.iterrows():
        donor_gt = row1['donor_gt']
        donor_querry = row1['donor_query']

        set1 =Genome_PiHAT[Genome_PiHAT['IID1'].str.contains(donor_querry)]
        set2 =Genome_PiHAT[Genome_PiHAT['IID2'].str.contains(donor_querry)]
        combo=pd.concat([set1,set2])

        set2_1 =combo[combo['IID1'].str.contains(donor_gt)]
        set2_2 =combo[combo['IID2'].str.contains(donor_gt)]
        combo2=pd.concat([set2_1,set2_2])

        try:
            GT_Match_Table.loc[i,f"PiHat: {Condition_Column}"]=combo2['PI_HAT'].values[0]
        except:
            _=''

    GT_Match_Table.to_csv(f"PiHAT_Stats_File_{Name}.csv",index=False,sep=',')


if __name__ == '__main__':
    main()
