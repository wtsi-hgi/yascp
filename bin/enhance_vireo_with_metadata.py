#!/usr/bin/env python3

__date__ = '2021-11-04'
__version__ = '0.0.1'

import argparse
import pandas as pd

parser = argparse.ArgumentParser(
    description="""
        Read an H5 10x-formatted object and write matrix files similar to
        10x output. This script requires pandas >1.0.1
        """
)

parser.add_argument(
    '-v', '--version',
    action='version',
    version='%(prog)s {version}'.format(version=__version__)
)

parser.add_argument(
    '-vireo_data', '--vireo_data',
    action='store',
    dest='vireo_data',
    default='vireo_exp__donor_n_cells.tsv',
)

parser.add_argument(
    '-em', '--Extra_Metadata_Donors',
    action='store',
    dest='Extra_Metadata_Donors',
    default='',
)

options = parser.parse_args()

pre_path=''
extra_metadata = pd.read_csv(options.Extra_Metadata_Donors,sep='\t')
vireo_data = pd.read_csv(options.vireo_data,sep='\t')
vireo_data = vireo_data.set_index('experiment_id')
for ix1 in vireo_data.index:
    tranche = ix1.split('__')[0]
    donor = ix1.split('__')[1]
    Assignemnts = pd.read_csv(f"GT_replace_{tranche}_assignments.tsv",sep='\t')
    Sample_replacement_data = Assignemnts[Assignemnts['donor_query']==donor]
    if (len(Sample_replacement_data)>0):
        if Sample_replacement_data['Match Expected'].values[0]:
            # print('yes')
            Sample_name = Sample_replacement_data['donor_gt'].values[0]
            split1 = Sample_name.split('___')
            if len(split1)==1:
                Sample_name =split1[0]
            else:
                Sample_name =split1[1]
            data_extra = extra_metadata[extra_metadata['experiment_id'] == f"{tranche}__{Sample_name}"]
            for col1 in data_extra.columns:
                if(col1!='experiment_id'):
                    try:
                        _=vireo_data.loc[f"{tranche}__{donor}"][col1]
                    except:
                        vireo_data[col1]='NA'
                    vireo_data.loc[f"{tranche}__{donor}",col1]=data_extra[col1].values[0]
# print(vireo_data)   
vireo_data.to_csv('replaced_vireo_exp__donor_n_cells_out.tsv',sep='\t')       

print('done')