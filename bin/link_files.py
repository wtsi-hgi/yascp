#!/usr/bin/env python3
import pandas as pd
import os
import glob
import argparse


parser = argparse.ArgumentParser(
    description="""
        Takes the required and optional inputs
        """
)


parser.add_argument(
    '-i', '--input_table',
    action='store',
    dest='input_table',
    required=True,
    help='input_table.'
)
options = parser.parse_args()
bridge_file = options.input_table

Data = pd.read_csv(bridge_file,sep='\t')
# Here we make the files to be in the right format before generating chanels.
Dataset3 = {}
c1 = 0
for i,row1 in Data.iterrows():
    print(i)
    c1 +=1
    # make a directory
    data_10x_format = row1.data_path_10x_format
    os.listdir(data_10x_format)
    outdir='input__'+row1.experiment_id
    if '/multi/count' in data_10x_format:
        data_10x_format = data_10x_format.split('/multi/count')[0]
    if os.path.isdir(data_10x_format+'/multi/count') or '/multi/count' in data_10x_format:
        # Cellranger 7
        os.mkdir(outdir)
        # Any vdj files?
        for vdj1 in glob.glob(data_10x_format+'/*/vdj*'):
            vdj1
            os.system(f"cd {outdir} && ln -s {vdj1} ./")
        # Metrics file:
        
        Metrics_file = glob.glob(data_10x_format+'/*/*/metrics_summary.csv')[0]
        Dataset = pd.read_csv(Metrics_file)
        Dataset1 = Dataset[Dataset['Library Type']=='Gene Expression']
        Dataset1 = Dataset1.set_index('Metric Name')
        Dataset1 = Dataset1['Metric Value']
        Dataset2 = pd.DataFrame(Dataset1).T
        Dataset2.to_csv(f'{outdir}/metrics_summary.csv',index=False)

        # ANalysis folder
        analysis_folder = glob.glob(data_10x_format+'/*/*/*/analysis')[0]
        os.symlink(analysis_folder, f"{outdir}/analysis")
        
        # filtered matrix
        sample_filtered_feature_bc_matrix = glob.glob(data_10x_format+'/*/*/count/sample_filtered_feature_bc_matrix')[0]
        os.symlink(f"{sample_filtered_feature_bc_matrix}", f"{outdir}/filtered_feature_bc_matrix")
        
        # raw_feature_bc_matrix
        raw_feature_bc_matrix =  glob.glob(data_10x_format+'/*/*count/raw_feature_bc_matrix')[0]
        os.symlink(raw_feature_bc_matrix, f"{outdir}/raw_feature_bc_matrix")
        
        #  sample_filtered_feature_bc_matrix
        sample_filtered_feature_bc_matrix =  glob.glob(data_10x_format+'/*/*/count/sample_filtered_feature_bc_matrix.h5')[0]
        os.symlink(sample_filtered_feature_bc_matrix, f"{outdir}/filtered_feature_bc_matrix.h5")
        
        # sample_molecule_info.h5
        sample_molecule_info = glob.glob(data_10x_format+'/*/*/count/sample_molecule_info.h5')[0]
        os.symlink(sample_molecule_info, f"{outdir}/molecule_info.h5")
        
        # sample_alignments.bam
        possorted_genome_bam = glob.glob(data_10x_format+'/*/*/count/sample_alignments.bam')[0]
        os.symlink(possorted_genome_bam, f"{outdir}/possorted_genome_bam.bam")
        
        # sample_alignments.bam.bai
        sample_alignments_bai = glob.glob(data_10x_format+'/*/*/count/sample_alignments.bam.bai')[0]
        os.symlink(sample_alignments_bai, f"{outdir}/possorted_genome_bam.bam.bai")
        
        # raw_feature_bc_matrix.h5
        raw_feature_bc_matrix = glob.glob(data_10x_format+'/*/count/raw_feature_bc_matrix.h5')[0]
        os.symlink(raw_feature_bc_matrix, f"{outdir}/raw_feature_bc_matrix.h5")
        
        # web_summary
        web_summary = glob.glob(data_10x_format+'/*/*/web_summary.html')[0]
        os.symlink(web_summary, f"{outdir}/web_summary.html")
        
    else:
        # Cellranger 6
        os.symlink(data_10x_format,outdir)
        
    Dataset3[c1]={'experiment_id':row1.experiment_id,'n_pooled':row1.n_pooled,'donor_vcf_ids':row1.n_pooled,'data_path_10x_format':os.getcwd()+'/'+outdir}
os.getcwd()
Yascp_input = pd.DataFrame(Dataset3).T  
Yascp_input.to_csv(f'input_file_corectly_formatted.tsv',sep='\t',index=False)

print('Done')              
