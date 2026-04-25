#!/usr/bin/env python

__date__ = '2022-06-30'
__version__ = '0.0.1'

import os
import pandas as pd
results_directoty = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/ELGH_26thMay_2022'
pool_id = results_directoty.split('/')[-1]
stream = os.popen(f'cd {results_directoty}/yascp && git log | head -n 1')
comit_tag = stream.read().replace('\n','')

stream = os.popen(f'cd {results_directoty}/yascp && git remote show origin | head -n2')
pipeline = stream.read().split('\n')[1].split(': ')[1]

stream = os.popen(f'ls -lt {results_directoty}/results/handover | head -n 2')
date_of_generation = stream.read()
dta = date_of_generation.split('\n')[1].split(' ')
date = '-'.join(dta[5:9])

Dataframe = pd.DataFrame({'pool_id':[pool_id],'pipeline':[pipeline],'comit_tag':[comit_tag],'date of run completion':[date],'results_directoty':[results_directoty]})
Dataframe.to_csv('pipeline_run_description.csv',sep='\t',index=False)