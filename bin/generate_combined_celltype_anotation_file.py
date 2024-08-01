#!/usr/bin/env python

__date__ = '2021-11-15'
__version__ = '0.0.1'

import argparse
from distutils.version import LooseVersion
import os
os.environ['NUMBA_CACHE_DIR']='/tmp'
os.environ['MPLCONFIGDIR']='/tmp'
import pandas as pd
import scanpy


def combine_reports(all_alternitive,mode):
    all_indexes_full=set({})
    for d1 in all_alternitive:
        if d1 in ('fake_file.fq', 'fake_file1.fq', 'fake_file2.fq'):
            Dataset = pd.DataFrame()
        else:
            Dataset = pd.read_csv(d1,sep='\t',index_col=0)
            if(len(Dataset.columns)==0):
                Dataset = pd.read_csv(d1,sep=',',index_col=0)
        Dataset=Dataset.add_prefix(mode)
        all_indexes = set(Dataset.index)
        all_indexes_full = all_indexes_full.union(all_indexes)
    Data_All_alt=pd.DataFrame(index=list(set(all_indexes_full)))    
    for d1 in all_alternitive:
        if d1 in ('fake_file.fq', 'fake_file1.fq', 'fake_file2.fq'):
            Dataset = pd.DataFrame()
        else:
            Dataset = pd.read_csv(d1,sep='\t',index_col=0)
            if(len(Dataset.columns)==0):
                Dataset = pd.read_csv(d1,sep=',',index_col=0)
        Dataset=Dataset.add_prefix(mode)
        for col1 in Dataset.columns:
            try:
                _ = Data_All_alt[col1]
            except:
                Data_All_alt[col1]=''
            Data_All_alt.loc[Dataset.index,col1] = Dataset[col1]
    return Data_All_alt

def main():
    """Run CLI."""
    parser = argparse.ArgumentParser(
        description="""
            Merges all the Celltype information in one file for Azimuth and Celltypist.
            """
    )

    parser.add_argument(
        '-v', '--version',
        action='version',
        version='%(prog)s {version}'.format(version=__version__)
    )

    parser.add_argument(
        '-az', '--all_azimuth_files',
        action='store',
        dest='all_azimuth_files',
        required=True,
        help='List of csv-delimited files of celltypes assigned by azimuth.'
    )
    
    parser.add_argument(
        '-af', '--all_other_paths',
        action='store',
        dest='all_alternitive',
        required=False,
        default=None,
        help='List of csv-delimited files of celltypes assigned by azimuth.'
    )

    parser.add_argument(
        '-ad', '--adata',
        action='store',
        dest='andata',
        required=True,
        help='Input adata to add labels to'
    )

    parser.add_argument(
        '-ct', '--all_celltypist_files',
        action='store',
        dest='all_celltypist_files',
        required=True,
        help='String of labels for each reduced_dims_tsv file. List should be\
            split by "::".'
    )

    options = parser.parse_args()
    
    
    Data_All=pd.DataFrame()
    
    # Read azimuth files from a TSV file using pandas
    azimuth_df = pd.read_csv(options.all_azimuth_files, header=None, names=['file_path'])
    azimuth_files = azimuth_df['file_path'].tolist()
    Data_All_Azimuth = combine_reports(azimuth_files,'Azimuth:')
    
    # Read celltypist files from a TSV file using pandas
    celltypist_df = pd.read_csv(options.all_celltypist_files, header=None, names=['file_path'])
    celltypist_files = celltypist_df['file_path'].tolist()

    celltypist_files2 = pd.DataFrame(celltypist_files,columns=['col1'])
    celltypist_files3 =list(celltypist_files2[~celltypist_files2['col1'].str.contains('input')]['col1'])
    Data_All_celltypist = combine_reports(celltypist_files3,'Celltypist:')
    
    if (options.all_alternitive):
        all_alternitive_df = pd.read_csv(options.all_alternitive, header=None, names=['file_path'])
        all_alternitive = all_alternitive_df['file_path'].tolist()
        Data_All_alt = combine_reports(all_alternitive,'')
    else:
        Data_All_alt=pd.DataFrame()

    Data_All = pd.concat([Data_All,Data_All_Azimuth,Data_All_celltypist,Data_All_alt],axis=1)
    
    Donor_Exp = Data_All.index.map(lambda x: '-'.join(x.split('-')[2:]))
    Donor = Donor_Exp.str.split('__').str[-1]
    Exp = Donor_Exp.str.split('__').str[0]

    Data_All['Donor'] =Donor
    Data_All['Exp'] =Exp
    Data_All.to_csv('All_Celltype_Assignments.tsv',sep='\t')

    adatas_df = pd.read_csv(options.andata, header=None, names=['file_path'])
    adatas = adatas_df['file_path'].tolist()
    adatasets = []
    # adatasets2 = adatasets[:2]
    adatasets__experiment_ids = []

    for ad1 in adatas:
        adata1 = scanpy.read_h5ad(ad1)
        if adata1.n_obs > 0:
            adatasets.append(adata1)

    ad = adatasets[0].concatenate(*adatasets[1:])
    if(len(adatasets)>1):
        # in this case the concentration adds a -1 -2 -3 to index that has to be removed.
        all_index = pd.DataFrame(ad.obs.index,columns=['col'])
        all_indexes = all_index['col'].str.split('-')
        all_together = all_indexes.apply(lambda x: '-'.join(x[:-1]))
        ad.obs.set_index(all_together, inplace=True)

        
    # ad2 = adatasets2[0].concatenate(*adatasets2[1:])
    # ad = scanpy.read(adata)
    ad.obs = ad.obs.merge(Data_All, left_index=True, right_index=True, how='left')
    # set(ad.obs.index)-set(Data_All.index)
    donor_celltype_report={}
    tranche_exp_report={}
    for id1 in set(Data_All['Exp']):
        print(id1)
        Exp_Data = Data_All[Data_All['Exp']==id1]
        
        dict_tranche_cells = {}
        for donor in set(Exp_Data['Donor']):
            dict_donor_cells = {}
            for col in Exp_Data.columns:
                if not 'score' in col and not 'probability' in col and not 'majority_voting' in col and not 'over_clustering' in col and not 'conf_score' in col and not 'n_counts' in col and not 'ver_clustering' in col and col !='Exp' and col !='Donor' and col !='n_counts':
                    print(col)
                    # col='Celltypist:over_clustering'
                    # col='Azimuth:predicted.celltype.l2'
                    counts = Exp_Data[col].value_counts()
                    counts.index = counts.index.astype(str) + ' - ' + col
                    dict_tranche_cells.update(counts.to_dict())
                    
                    # print(donor)
                    donor_data=Exp_Data[Exp_Data['Donor']==donor]
                    donor_counts = donor_data[col].value_counts()
                    donor_counts.index = donor_counts.index.astype(str) + ' - '+col
                    dict_donor_cells.update(donor_counts.to_dict())
                # check all the available celltypes here
                # and count the numbers
            donor_celltype_report[f'{id1} {donor}']=dict_donor_cells
        tranche_exp_report[id1]=dict_tranche_cells
        

    pd_tranche_exp_report = pd.DataFrame(tranche_exp_report).T
    pd_tranche_exp_report =pd_tranche_exp_report.fillna(0)
    # Generate Donor Report of counts
    pd_donor_celltype_report = pd.DataFrame(donor_celltype_report).T
    pd_donor_celltype_report =pd_donor_celltype_report.fillna(0)

    pd_donor_celltype_report.to_csv('donor_celltype_report.tsv',sep='\t')
    pd_tranche_exp_report.to_csv('tranche_celltype_report.tsv',sep='\t')
    
    # Generate Experiment report of counts
    ad.write(
        'adata.h5ad',
        compression='gzip',
        compression_opts=9  # takes ages, but we want a small file for system
    )

if __name__ == '__main__':
    main()
