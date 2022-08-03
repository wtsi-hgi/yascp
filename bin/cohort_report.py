#!/usr/bin/env python3
import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(
    description="""
        Combines all the Raw cellranger metadata to be passed in the h5ad
        """
)

parser.add_argument(
    '-d',
    action='store',
    dest='path',
    required=False,
    default='./results',
    help='Path')

options = parser.parse_args()
path = options.path

# Split Donor Report by cohort
project_name = os.listdir(f"{path}/handover/Summary_plots")[0]
GT_MATCH = pd.read_csv(f"{path}/deconvolution/vireo_gt_fix/assignments_all_pools.tsv",sep='\t')
Donor_Report = pd.read_csv(f"{path}/handover/Donor_Quantification_summary/Donor_Report.tsv",sep='\t')
Tranch_Report = pd.read_csv(f"{path}/handover/Donor_Quantification_summary/Tranche_Report.tsv",sep='\t')
Extra_Metadata_Donors = pd.read_csv(f"{path}/handover/Summary_plots/{project_name}/Summary/Extra_Metadata_Donors.tsv",sep='\t')
# qc/Cardinal_45327_Jul_18_2022/work/62/aacd091ea7711fae264e884e122230/Summary_plots/Cardinal_45327_Jul_18_2022/Summary
# t2 = GT_MATCH.loc[GT_MATCH['donor_gt'].str.contains('celline')]
# GT_MATCH.loc[GT_MATCH['donor_gt'].str.contains('celline'),'Match Expected']=True
GT_MATCH_CONFIDENT = GT_MATCH[GT_MATCH['Match Expected'] == True]

def Generate_Report(GT_MATCH_CONFIDENT,pan):

    Total_Report = pd.DataFrame()               
    for i,row1 in GT_MATCH_CONFIDENT.iterrows():
        print(i)
        # print(row1)
        print(f"{row1['pool']}_{row1['donor_query']}")
        
        Matched_Donor_report = ''
        Matched_Donor_report = Donor_Report[Donor_Report['Pool_ID.Donor_Id']==f"{row1['pool']}_{row1['donor_query']}"]
        
        Tranch_stats = Tranch_Report[Tranch_Report['Pool id'] ==row1['pool']]
        for i,row12 in Matched_Donor_report.iterrows():
            print(i)
            print(row12)
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
        cohort = row1['final_panel']
        try:
            live_cell_count = Extra_Metadata_Donor1['live_cell_count'].values[0]
            viability =Extra_Metadata_Donor1['viability'].values[0]
            SITE = Extra_Metadata_Donor1['SITE'].values[0]
            Amount = Extra_Metadata_Donor1['Amount'].values[0]
            RECIEVED = Extra_Metadata_Donor1['RECIEVED'].values[0]
        except:
            live_cell_count = 'NA'
            viability ='NA'
            SITE = 'NA'
            Amount = 'NA'
            RECIEVED = 'NA'

        Matched_Donor_report.insert(8, 'lab_live_cell_count',live_cell_count)
        Matched_Donor_report.insert(8, 'viability', viability)
        Matched_Donor_report.insert(8, 'cohort', cohort)
        Matched_Donor_report.insert(8, 'site', SITE)
        Matched_Donor_report.insert(8, 'amount recieved', Amount)
        Matched_Donor_report['Date sample received'] =  RECIEVED
        Matched_Donor_report['Donor id']=row1['donor_gt original'].split('_')[0]
        Total_Report = pd.concat([Total_Report,Matched_Donor_report])

    return Total_Report


for confident_panel in set(GT_MATCH['final_panel']):
    if (confident_panel!='NONE' and confident_panel!='GT_cell_lines'):
        print(confident_panel)
        pan=confident_panel.replace('GT_','')
        print(pan)
        GR_PANEL = GT_MATCH[GT_MATCH['final_panel'] == confident_panel]
        GT_CELLINE = GT_MATCH[GT_MATCH['final_panel'] == 'GT_cell_lines']
        GT_PANEL_CELLINE = pd.concat([GR_PANEL,GT_CELLINE])
        Total_Report = Generate_Report(GT_PANEL_CELLINE,pan)
        Expected_Samples = Extra_Metadata_Donors[Extra_Metadata_Donors.cohort == confident_panel]
        Missing_Samples = set(Expected_Samples.donor)-set(Total_Report['Vacutainer ID'])
        try:
            os.mkdir(f'{path}/handover/Summary_plots/{project_name}/Summary/{pan}_REPORT')
        except:
            print('Dir exists')
        if (len(Missing_Samples)>0):
            Missing = Extra_Metadata_Donors.loc[Missing_Samples]['experiment_id']
            Missing.to_csv(f'{path}/handover/Summary_plots/{project_name}/Summary/{pan}_REPORT/{project_name}_Missing_UKB_Donors.tsv',sep='\t')
        Total_Report.to_csv(f'{path}/handover/Summary_plots/{project_name}/Summary/{pan}_REPORT/{project_name}_UKBB_Report.tsv',sep='\t',index=False)

print('Done')