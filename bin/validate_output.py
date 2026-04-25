#!/usr/bin/env python3

import pandas as pd
import glob
from functools import reduce
import sys

def vireo_folders_processing ():
    #collect number of total barcodes, failed barcodes and donors per pool
    filepath_list=glob.glob('vireo_*/summary.tsv')
    pool_list=[filepath.lstrip('vireo_').rstrip('/summary.tsv') for filepath in filepath_list]
    vireo_res=[vireo_summary_processing(filepath) for filepath in filepath_list]
    donor_counts_list=[res[0] for res in vireo_res]
    total_cell_counts=[res[1] for res in vireo_res]
    unasigned_and_duplets=[res[2] for res in vireo_res]
    df = pd.DataFrame({
    'Pool': pool_list,
    'donor_counts': donor_counts_list,
    'total_cell_counts': total_cell_counts,
    'unasigned_and_duplets': unasigned_and_duplets
    })
    return df

def vireo_summary_processing (input_file):
    #process vireo summary.tsv for each pool
    df=pd.read_csv(input_file, sep="\t")
    donor_counts=len(df)-2
    total_cell_counts=df['Freq'].sum()
    df.set_index('Var1', inplace=True)
    unasigned_and_duplets=df.at['doublet','Freq']+df.at['unassigned','Freq']

    return [donor_counts, total_cell_counts, unasigned_and_duplets]

def read_counts (file_path):
    with open(file_path, "r") as f:
        number = int(f.read().strip().split()[0])

    return number

def collect_counts_per_pool (suffix, count_name):
    #collect number of barcodes per pool
    file_list=glob.glob('*'+suffix)
    pool_dict={}
    for file_name in file_list:
        pool=file_name.replace(suffix, '')
        pool_dict[pool]=read_counts(file_name)
    pool_df = pd.DataFrame(list(pool_dict.items()), columns=["Pool", count_name])
    return pool_df

def process_pca(file_name, donor_counts_title, pca_barcode_counts_title):
    #process adata-pcs.tsv.gz to collect number of barcodes and donors per pool processed in  clustering and integration
    df=pd.read_csv(file_name, sep="\t", compression="gzip")
    df=df[['cell_barcode']]
    #pattern = r'^([ACGT]+-\d+)-(.+)__(donor\d+)$'
    pattern = r'^(.+)-([^-]+)__(donor\d+)$'
    df[['barcode', 'Pool', 'donor']] = df['cell_barcode'].str.extract(pattern)
    donor_counts = df.groupby('Pool')['donor'].nunique().reset_index(name=donor_counts_title)
    barcode_counts = df.groupby('Pool')['barcode'].nunique().reset_index(name=pca_barcode_counts_title)
    merged_df=pd.merge(donor_counts, barcode_counts, on="Pool", how="outer")
    return merged_df

def main ():
    azimuth_run = sys.argv[1]
    celltypist_run = sys.argv[2]
    scpred_run = sys.argv[3]

    sorting_list1=['Pool', 'preprocessing']
    sorting_list2=[]
    sorting_list3=[]
    sorting_list4=[]
    #preprocess
    preprocessing_by_pool_df=collect_counts_per_pool('_barcodes.counts.txt', 'preprocessing')
    df=preprocessing_by_pool_df
    if glob.glob("*_barcodes_filtered.counts.txt"):
        preprocessing_by_pool_filtered_df=collect_counts_per_pool('_barcodes_filtered.counts.txt', 'preprocessing_filtered')
        df=pd.merge(preprocessing_by_pool_df, preprocessing_by_pool_filtered_df, on="Pool", how="outer")
        sorting_list1.append('preprocessing_filtered')
        sorting_list3=['preprocessing_filtered']
    else:
        prpocess_cellbender=collect_counts_per_pool('_cellbender.counts.txt', 'preprocessing_filtered')
        df=pd.merge(preprocessing_by_pool_df, prpocess_cellbender, on="Pool", how="outer")
        sorting_list1.append('preprocessing_filtered')
        sorting_list3=['preprocessing_filtered']

    #dedup
    if glob.glob("*_doublet.counts.txt"):
        doublets_by_pool_df=collect_counts_per_pool('_doublet.counts.txt', 'doublet_detection')
        df=pd.merge(df, doublets_by_pool_df, on="Pool", how="outer")
        sorting_list1.append('doublet_detection')
        sorting_list3.append('doublet_detection')

    #celltype
    if glob.glob("cells_by_pool.counts.txt"):
        cells_by_pool_df= pd.read_csv('cells_by_pool.counts.txt', sep="\t", header=None, names=["Pool","Azimuth","celltypist","All_alt"])
        df=pd.merge(df, cells_by_pool_df, on="Pool", how="outer")
        if azimuth_run == "true":
            sorting_list1.extend(['Azimuth'])
            sorting_list3.extend(['Azimuth'])
        if celltypist_run == "true":
            sorting_list1.extend(['celltypist'])
            sorting_list3.extend(['celltypist'])
        if scpred_run == "true":
            sorting_list1.extend(['All_alt'])
            sorting_list3.extend(['All_alt'])

    #deconvolution
    if glob.glob("*_cellSNP.counts.txt"):
        cellSNP_by_pool_df=collect_counts_per_pool('_cellSNP.counts.txt', 'cellSNP')
        df=pd.merge(df, cellSNP_by_pool_df, on="Pool", how="outer")
        sorting_list1.append('cellSNP')
        sorting_list3.append('cellSNP')
        vireo_df=vireo_folders_processing()
        df=pd.merge(df, vireo_df, on="Pool", how="outer")
        sorting_list1.extend(['total_cell_counts', 'unasigned_and_duplets'])
        sorting_list3.extend(['total_cell_counts'])
        sorting_list2.append('donor_counts')
        genotypes_by_pool_df=collect_counts_per_pool('_infered_genotypes.counts.txt', 'infered_genotypes')
        df=pd.merge(df, genotypes_by_pool_df, on="Pool", how="outer")
        sorting_list2.append('infered_genotypes')

    #introduce genotype matcher

    #clustering and integration
    if glob.glob("adata-pcs.tsv.gz"):
        pca_df=process_pca('adata-pcs.tsv.gz', 'pca_donor_counts', 'pca_barcode_counts')
        df=pd.merge(df, pca_df, on="Pool", how="outer")
        sorting_list1.append('pca_barcode_counts')
        #sorting_list3.append('pca_barcode_counts')
        sorting_list2.append('pca_donor_counts')
        sorting_list4.append('pca_barcode_counts')

    if glob.glob("*-reduced_dims.tsv.gz"):
        harmony_files = glob.glob("*-reduced_dims.tsv.gz")
        harmony_df=process_pca(harmony_files[0], 'harmony_donor_counts', 'harmony_barcode_counts')
        df=pd.merge(df, harmony_df, on="Pool", how="outer")
        sorting_list1.append('harmony_barcode_counts')
        #sorting_list3.append('harmony_barcode_counts')
        sorting_list2.append('harmony_donor_counts')
        sorting_list4.append('harmony_barcode_counts')

    if glob.glob("reduced_dims-*.tsv.gz"):
        bbknn_files = glob.glob("reduced_dims-*.tsv.gz")
        bbknn_df=process_pca(bbknn_files[0], 'bbknn_donor_counts', 'bbknn_barcode_counts')
        df=pd.merge(df, bbknn_df, on="Pool", how="outer")
        sorting_list1.append('bbknn_barcode_counts')
       #sorting_list3.append('bbknn_barcode_counts')
        sorting_list2.append('bbknn_donor_counts')
        sorting_list4.append('bbknn_barcode_counts')

    #introduce bbkn
    #######################################

    if sorting_list2:
        sorting_list1.extend(sorting_list2)
    df = df.reindex(sorting_list1, axis=1)
    df.to_csv("validation_data_combined.tsv", sep="\t", index=False)


    if 'total_cell_counts' in sorting_list1:
        df['passed_cell_counts']=df['total_cell_counts']-df['unasigned_and_duplets']
    
    if (len(sorting_list3)>1) and ('preprocessing_filtered' in sorting_list3):
        df['barcode_counts_consistency'] = df[sorting_list3].eq(df['preprocessing_filtered'], axis=0).all(axis=1)

    if len(sorting_list2)>1:
        df['donor_counts_consistency'] = df[sorting_list2].eq(df['donor_counts'], axis=0).all(axis=1)

    if ('total_cell_counts' in sorting_list1) and len(sorting_list4) >0:
        sorting_list4.append('passed_cell_counts')
        df['passed_barcode_counts_consistency'] =  df[sorting_list4].eq(df['passed_cell_counts'], axis=0).all(axis=1)


    consistency_cols=['barcode_counts_consistency','passed_barcode_counts_consistency', 'donor_counts_consistency']
    if any(col in df.columns for col in consistency_cols):
        with open ("results_validation_report.txt", "w") as output_res:
            if 'passed_barcode_counts_consistency' in df.columns:
                if not df['barcode_counts_consistency'].all() or not df['passed_barcode_counts_consistency'].all():
                    output_res.write("INCONSISTENCY IN BARCODE NUMBERS!\n")
                else:
                    output_res.write("Barcode numbers are consistent\n")
            elif 'barcode_counts_consistency' in df.columns:
                if not df['barcode_counts_consistency'].all():
                    output_res.write("INCONSISTENCY IN BARCODE NUMBERS!\n")
                else:
                    output_res.write("Barcode numbers are consistent\n")
            if 'donor_counts_consistency' in df.columns:
                if not df['donor_counts_consistency'].all():
                    output_res.write("INCONSISTENCY IN DONOR NUMBERS!\n")
                else:
                    output_res.write("Donor numbers are consistent\n")

if __name__ == '__main__':
    main()