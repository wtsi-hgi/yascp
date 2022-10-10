import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Replace Donor Ids for GT match")
parser.add_argument("-e", "--exp_id", metavar="exp_id", dest = "exp_id",
                    help="", type=str,
                    default=1)
parser.add_argument("-res", "--res", metavar="res", dest = "res",
                    help="", type=str,
                    default='results')

parser.add_argument("-met", "--met_folder", metavar="met_folder", dest = "met_folder",
                    help="", type=str,
                    default='results')

args = parser.parse_args()
exp_id = args.exp_id
Report_Donor = pd.read_csv(f'results/handover/Donor_Quantification_summary/{exp_id}_Donor_Report.tsv',sep='\t')
Tranche_report = pd.read_csv(f'results/handover/Donor_Quantification_summary/{exp_id}_Tranche_Report.tsv',sep='\t')

Extra_Metadata = pd.read_csv(f'{args.met_folder}',sep='\t')
Extra_Metadata = Extra_Metadata.set_index('experiment_id')
Report_Donor['Pool ID']

for unique_pool_id in set(Extra_Metadata.index):
    Report_Donor.loc[Report_Donor['Pool ID']==unique_pool_id,'Date of sample sequencing']=Extra_Metadata.loc[unique_pool_id,'last_updated']
    Tranche_report.loc[Tranche_report['Pool id']==unique_pool_id,'Date of sample sequencing']=Extra_Metadata.loc[unique_pool_id,'last_updated']
    Tranche_report.loc[Tranche_report['Pool id']==unique_pool_id,'UKB donors expected in pool']=Extra_Metadata.loc[unique_pool_id,'nr_ukbb_samples']
    Tranche_report.loc[Tranche_report['Pool id']==unique_pool_id,'ELGH donors expected in the pool']=Extra_Metadata.loc[unique_pool_id,'nr_elgh_samples']
    Tranche_report.loc[Tranche_report['Pool id']==unique_pool_id,'Spikeins expected in the pool']=Extra_Metadata.loc[unique_pool_id,'nr_spikeins']

Tranche_report.to_csv(f'results/handover/Donor_Quantification_summary/{exp_id}_Tranche_Report.tsv',sep='\t',index=False)
Report_Donor.to_csv(f'results/handover/Donor_Quantification_summary/{exp_id}_Donor_Report.tsv',sep='\t',index=False)
# now load the updated metadata

print('Done')