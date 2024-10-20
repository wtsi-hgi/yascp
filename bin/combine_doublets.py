#!/usr/bin/env python

import pandas
import argparse
import pandas as pd

# defined command line options
# this also generates --help and error handling
CLI=argparse.ArgumentParser()
CLI.add_argument(
  "--list",  # name on the CLI - drop the `--` for positional/required parameters
  nargs="*",  # 0 or more values expected => creates a list
  type=str,
  default=None,  # default if nothing is provided
)

args = CLI.parse_args()
all_files = args.list
all_combo = pd.DataFrame()
for f1 in all_files:
    print(f1)
    df1 = pd.read_csv(f1, index_col=0, sep='\t')
    all_combo = pd.concat([all_combo, df1], axis=1)
all_combo.index.name = 'barcodes'
try:
  all_combo['Scrublet_DropletType'] = all_combo['scrublet__predicted_multiplet']
  all_combo.loc[all_combo['Scrublet_DropletType'] == True, 'Scrublet_DropletType'] = 'doublet'
  all_combo.loc[all_combo['Scrublet_DropletType'] == False, 'Scrublet_DropletType'] = 'singlet'
except:
  print('Scrublet wast exacuted')
all_combo.to_csv('all_doublet_results_combined.tsv',sep='\t') 
print('Done')