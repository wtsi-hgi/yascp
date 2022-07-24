#!/usr/bin/env python3

import pandas as pd
import os
# Split Donor Report by cohort
project_name = os.listdir("./Summary_plots")[0]
GT_MATCH = pd.read_csv("./results/deconvolution/vireo_gt_fix/assignments_all_pools.tsv",sep='\t')
Donot_Report = pd.read_csv("./results/handover/Donor_Quantification_summary/Donor_Report.tsv",sep='\t')
Tranch_Report = pd.read_csv("./results/handover/Donor_Quantification_summary/Tranche_Report.tsv",sep='\t')
Extra_Metadata_Donors = pd.read_csv(f"./Summary_plots/{project_name}/Summary/Extra_Metadata_Donors.tsv",sep='\t')
# qc/Cardinal_45327_Jul_18_2022/work/62/aacd091ea7711fae264e884e122230/Summary_plots/Cardinal_45327_Jul_18_2022/Summary
# t2 = GT_MATCH.loc[GT_MATCH['donor_gt'].str.contains('celline')]
GT_MATCH.loc[GT_MATCH['donor_gt'].str.contains('celline'),'Match Expected']=True
GT_MATCH_CONFIDENT = GT_MATCH[GT_MATCH['Match Expected'] == True]

def Generate_Report(GT_MATCH_CONFIDENT):

    Total_Report = pd.DataFrame()               
    for i,row1 in GT_MATCH_CONFIDENT.iterrows():
        print(i)
        print(row1)
        print(f"{row1['pool']}_{row1['donor_query']}")
        
        Matched_Donor_report = ''
        Matched_Donor_report = Donot_Report[Donot_Report['Donor id_Experiment ID']==f"{row1['donor_query']}_{row1['pool']}"]
        
        Tranch_stats = Tranch_Report[Tranch_Report['Pool id'] ==row1['pool']]
        # for i,row12 in Tranch_stats.iterrows():
        #     print(i)
        #     print(row12)
        if ('celline' in row1['donor_gt']):
            print('celline')
            gt_match = row1['donor_gt'].split('___')[1].replace('celline_','')
        else:
            gt_match = row1['donor_gt']
        
        Extra_Metadata_Donor1 = Extra_Metadata_Donors[Extra_Metadata_Donors['experiment_id'].str.contains(f"{row1['pool']}__{gt_match}")]
        
        Matched_Donor_report.insert(2, "Vacutainer ID", gt_match)
        Matched_Donor_report['Experiment ID'] = project_name
        Matched_Donor_report['Chromium channel number'] = Tranch_stats['Chromium channel number'].values[0]
        Matched_Donor_report['Date of sample sequencing'] = Tranch_stats['Date of sample sequencing'].values[0]
        Matched_Donor_report.insert(8, 'lab_live_cell_count', Extra_Metadata_Donor1['live_cell_count'].values[0])
        Matched_Donor_report.insert(8, 'viability', Extra_Metadata_Donor1['viability'].values[0])
        Matched_Donor_report.insert(8, 'cohort', Extra_Metadata_Donor1['cohort'].values[0])
        Matched_Donor_report.insert(8, 'site', Extra_Metadata_Donor1['SITE'].values[0])
        Matched_Donor_report.insert(8, 'amount recieved', Extra_Metadata_Donor1['Amount'].values[0])
        Matched_Donor_report['Date sample received'] =  Extra_Metadata_Donor1['RECIEVED'].values[0]
        Total_Report = pd.concat([Total_Report,Matched_Donor_report])
        # for i,row12 in t2.iterrows():
        #     print(i)
        #     print(row12)
        # for i,row12 in GT_MATCH.iterrows():
        #     print(i)
        #     print(row12)
    return Total_Report

GT_ELGH = GT_MATCH[GT_MATCH['final_panel'] == 'GT_ELGH']
GT_CELLINE = GT_MATCH[GT_MATCH['final_panel'] == 'GT_cell_lines']
GT_ELGH_CELLINE = pd.concat([GT_ELGH,GT_CELLINE])

GT_ELGH_SUC = GT_ELGH[GT_ELGH['Match Expected'] == True]

# Total_Report_ELGH = Generate_Report(GT_ELGH_CELLINE)

Total_Report = Generate_Report(GT_MATCH_CONFIDENT)
Total_Report_UKBB_Samples = Total_Report[Total_Report.cohort == 'Cardinal UKB']
Total_Report_Controls = Total_Report[Total_Report.cohort == 'Cardinal controls']
Total_Report2=pd.concat([Total_Report_Controls,Total_Report_UKBB_Samples])
UKBB_Samples = Extra_Metadata_Donors[Extra_Metadata_Donors.cohort == 'Cardinal UKB']
Missing_UKBB = set(UKBB_Samples.donor)-set(Total_Report2['Vacutainer ID'])
Extra_Metadata_Donors2 = Extra_Metadata_Donors.set_index('donor') 
Missing = Extra_Metadata_Donors2.loc[Missing_UKBB]['experiment_id']
os.mkdir(f'Summary_plots/{project_name}/Summary/UKBB_REPORT')
Missing.to_csv(f'Summary_plots/{project_name}/Summary/UKBB_REPORT/{project_name}_Missing_UKB_Donors.tsv')
Total_Report2.to_csv(f'Summary_plots/{project_name}/Summary/UKBB_REPORT/{project_name}_UKBB_Report.tsv',sep='\t',index=False)

# Most of the ELGH samples are not expected as they provide ungentyped donors


print('Done')