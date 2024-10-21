#!/usr/bin/env python
import glob
import pandas as pd
f2=glob.glob('./*.barcodes.tsv')

combo = pd.DataFrame()
for f in f2:
    don = f.split('/')[-1].split('.')[0]
    D2= pd.read_csv(f,sep='\t',header=None,names=['cell'])
    D2['donor_id']=don
    combo = pd.concat([combo,D2])
combo.to_csv('donor_ids.tsv',sep='\t',index=False)