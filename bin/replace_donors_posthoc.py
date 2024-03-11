#!/usr/bin/env python3
import pandas as pd
import argparse


parser = argparse.ArgumentParser(description="Replace Donor Ids for GT match")
parser.add_argument("-i", "--in_file", metavar="donor_id", dest = "in_file",
                    help="", type=str,
                    default=1)
parser.add_argument("-gp", "--genotype_phenotype_mapping", metavar="genotype_phenotype_mapping", dest = "genotype_phenotype_mapping",
                    help="", type=str,
                    default=None)

parser.add_argument("-if", "--input_file", metavar="input_file", dest = "input_file",nargs='+',
                    help="", type=str,
                    default=None)
                    
args = parser.parse_args()
if (args.genotype_phenotype_mapping):
    gp_ma = pd.read_csv(args.genotype_phenotype_mapping,sep='\t',index_col=0)
    gp_ma.index = gp_ma.index.astype(str)
else:
    gp_ma = pd.DataFrame()
input_id =args.in_file

try:
    GT_Assignments = pd.read_csv(f'{input_id}',sep=',')
    if len(GT_Assignments.columns)==1:
        GT_Assignments = pd.read_csv(f'{input_id}',sep='\t')
except:
    GT_Assignments = pd.read_csv(f'{input_id}',sep='\t')

for i,row1 in GT_Assignments.iterrows():
    # print(row1)
    
    replacement=row1['donor_gt original']

    expected = str(row1['Emergency_ids expected'])+','+str(row1['Good_ids expected'])
    expected= expected.split(',')
    D2 = pd.DataFrame(expected, columns=['col1'])
    # try:
    #     repl = str(gp_ma.loc[math.split('_')[0]][0])
    #     GT_Assignments.loc[i,'donor_gt']=repl
    #     if (len(DF2[DF2.col1.str.contains(repl)])>0):
    #         GT_Assignments.loc[i,'Match Expected']=True
    # except:
    #     try:
    #         repl = str(gp_ma.loc[math][0])
    #         GT_Assignments.loc[i,'donor_gt']=repl
    #         if (len(DF2[DF2.col1.str.contains(repl)])>0):
    #             GT_Assignments.loc[i,'Match Expected']=True
    #     except:
    #         print('no need to replace')
            
            
            
            #####
            
    try:
        # replacement = 'S2-046-01741'
        try:
            replacement = str(gp_ma.loc[replacement].iloc[0].values[0])
        except:
            replacement = str(gp_ma.loc[replacement].iloc[0])
        # print(replacement)
        # replacement = replacement.values[0]
        
        if len(D2[D2.col1.str.contains(replacement)])>0:
            GT_Assignments.loc[i,'Match Expected']='True'
            remove_from_set = D2[D2.col1.str.contains(replacement)]['col1'].values[0]
            All_Expected_set.append(remove_from_set)
    except:
        try:
            # replacement='15001608190388_204238910153_R09C02'
            
            try:
                replacements = gp_ma.loc[replacement.split('_')[0]]
            except:
                try:
                    replacements = gp_ma.loc[replacement.split('_')[1]]
                except:
                    replacements = gp_ma.loc[replacement.split('_')[2]]
            
            if len(replacements)>1:
                replacement = ''
                for rep1 in replacements.iloc[:,0]:
                    
                    if len(D2[D2.col1.str.contains(rep1)])>0:
                        GT_Assignments.loc[i,'Match Expected']='True'
                        replacement = rep1
                        break
                        remove_from_set = D2[D2.col1.str.contains(replacement)]['col1'].values[0]
                        All_Expected_set.append(remove_from_set)
                    else:
                        replacement=replacement+rep1+';'
            else:
                replacement = replacements.values[0]
        except:
            replacement = 'No_mapping___'+replacement
            
    GT_Assignments.loc[i,'donor_gt']=replacement
GT_Assignments.to_csv(input_id,sep='\t',index=False)