import pandas as pd
# Split Donor Report by cohort

GT_MATCH = pd.read_csv("/lustre/scratch123/hgi/projects/cardinal_analysis/qc/UKBB_ELGH_5th_July_2022/Summary_plots/UKBB_ELGH_5th_July_2022/Summary/assignments_all_pools.tsv",sep='\t')
Donot_Report = pd.read_csv("/lustre/scratch123/hgi/projects/cardinal_analysis/qc/UKBB_ELGH_5th_July_2022/Summary_plots/UKBB_ELGH_5th_July_2022/Summary/Donor_Report.tsv",sep='\t')
Tranch_Report = pd.read_csv("/lustre/scratch123/hgi/projects/cardinal_analysis/qc/UKBB_ELGH_5th_July_2022/Summary_plots/UKBB_ELGH_5th_July_2022/Summary/Tranche_Report_UKBB_ELGH_5th_July_2022.tsv",sep='\t')
Extra_Metadata_Donors = pd.read_csv("/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/fetch/UKBB_ELGH_5th_July_2022/results/yascp_inputs/Extra_Metadata_Donors.tsv",sep='\t')

Total_Report = pd.DataFrame()
GT_MATCH_CONFIDENT = GT_MATCH[GT_MATCH['Match Expected'] == True]
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
    Extra_Metadata_Donor1 = Extra_Metadata_Donors[Extra_Metadata_Donors['experiment_id']==f"{row1['pool']}__{row1['donor_gt']}"]
    
    Matched_Donor_report.insert(2, "Vacutainer ID", row1['donor_gt'])
    Matched_Donor_report['Experiment ID'] = Tranch_stats['Experiment id'].values[0]
    Matched_Donor_report['Chromium channel number'] = Tranch_stats['Chromium channel number'].values[0]
    Matched_Donor_report['Date of sample sequencing'] = Tranch_stats['Date of sample sequencing'].values[0]
    Matched_Donor_report.insert(8, 'lab_live_cell_count', Extra_Metadata_Donor1['live_cell_count'].values[0])
    Matched_Donor_report.insert(8, 'viability', Extra_Metadata_Donor1['viability'].values[0])
    Matched_Donor_report.insert(8, 'cohort', Extra_Metadata_Donor1['cohort'].values[0])
    Matched_Donor_report.insert(8, 'site', Extra_Metadata_Donor1['SITE'].values[0])
    Matched_Donor_report.insert(8, 'amount recieved', Extra_Metadata_Donor1['Amount'].values[0])
    Matched_Donor_report['Date sample received'] =  Extra_Metadata_Donor1['RECIEVED'].values[0]
    Total_Report = pd.concat([Total_Report,Matched_Donor_report])
    # for i,row12 in Matched_Donor_report.iterrows():
    #     print(i)
    #     print(row12)
    # for i,row12 in Extra_Metadata_Donor1.iterrows():
    #     print(i)
    #     print(row12)
# Missing_Vacutainers =
Total_Report.to_csv('/lustre/scratch123/hgi/projects/cardinal_analysis/qc/UKBB_ELGH_5th_July_2022/Summary_plots/UKBB_ELGH_5th_July_2022/Summary/UKBB_Repoty.tsv',sep='\t',index=False)
print('Done')