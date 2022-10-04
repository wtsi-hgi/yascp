#!/usr/bin/env python

__date__ = '2022-09-20'
__version__ = '0.0.1'

import argparse
import pandas as pd
import re

def main():
    print('main')
    parser = argparse.ArgumentParser(
        description="""
            Takes the required and optional inputs
            """
    )
    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-b', '--bridge',
        action='store',
        dest='bridge',
        required=False,
        default=None,
        help='Bridge file to use to convert samples to the mappings'
    )

    parser.add_argument(
        '-s', '--samples',
        action='store',
        dest='samples',
        required=True,
        help='Samples.'
    )

    parser.add_argument(
        '-vs', '--vcfsamples',
        action='store',
        dest='vcfsamples',
        required=True,
        help='vcfsamples.'
    )

    parser.add_argument(
        '-o', '--output',
        action='store',
        dest='output',
        default='output.csv',
        required=False,
        help='Output.'
    )

    options = parser.parse_args()
    bridge_file = options.bridge
    samples_string = options.samples
    vcfsamples = options.vcfsamples
    output_file = options.output

    # if bridge file is provided we will attempt to map it to the new entries, otherwise keep the same. 
    if bridge_file:
        print('yes')
        bridge = pd.read_csv(bridge_file,sep='\t')
    else:
        bridge = pd.DataFrame(columns=['oragene_id','s00046_id'])

    samples = samples_string.split(',')
    all_maped_samples = []
    for s1 in samples:
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
        all_maped_samples.append(replacement)

    # Determine the overlap between samples and available entries in the vcf file
    Data_vcfsamples = pd.read_csv(vcfsamples,header=None)
    Data_vcfsamples.columns=['samples']
    # Interset = set(all_maped_samples).intersection(set(Data_vcfsamples.iloc[:,0]))
    Interset=set()
    # we look fo partial patterns in dataset as ELGH samples is typically designed in this way.
    # example mapping - 15001506221453 - full id = 
    for m1 in all_maped_samples:
        d2 = Data_vcfsamples[Data_vcfsamples['samples'].str.contains(m1)]
        if len(d2)>0:
          i2 = d2.values[0][0]
          Interset.add(i2)
    Interset = pd.DataFrame(Interset,columns=['samples'])
    if len(Interset)>0:
        Interset=Interset.sort_values(by=['samples'])
        Interset.to_csv(output_file,index=False,header=False)
    print('Done')





if __name__ == '__main__':
    main()