#!/usr/bin/env python

import pandas as pd
import click
import argparse
import scanpy as sc

parser = argparse.ArgumentParser(description="Replace Donor Ids for GT match")
parser.add_argument("-s", "--scrublet", metavar="scrublet", dest = "scrublet",
                    help="", type=str,
                    default=None)
args = parser.parse_args()

# scrublet = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/carls_data/analysis_with_gt/work/4b/112b9aaf014e97a8dac1b9e52f210b/WC101549143019__doublet_results_combined.tsv"
if __name__ == '__main__':
    try:
        scrublet_data = pd.read_csv(args.scrublet,sep='\t')
    except:
        scrublet_data = sc.read_10x_mtx(args.scrublet)
        scrublet_data = pd.DataFrame(scrublet_data.obs.index,columns=['barcodes'])
        scrublet_data['donor_id']='donor'
    scrublet_data = scrublet_data.rename(columns={'barcodes':'cell','Scrublet_DropletType':'donor_id'})
    scrublet_data.loc[scrublet_data['donor_id']=='singlet','donor_id'] = 'donor'
    scrublet_data.to_csv('cell_belongings.tsv',sep='\t',index=False)
    print("Done")