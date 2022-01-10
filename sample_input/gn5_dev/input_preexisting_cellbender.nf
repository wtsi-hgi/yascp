params {

    run_gather_data = true
    
    input = 'existing_cellbender'
    cellbender_file='/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/cellbender/Deconvolution_Exp4/outputs/qc_cluster_input_files/file_paths_10x-cellbender_params__epochs_250__learnrt_1pt0Eneg7__zdim_100__zlayer_500__lowcount_10-FPR_0pt1.tsv' //if cellbender is run already then can skip this by selecting existing_cellbender and input 
    extra_metadata = ''
    skip_preprocessing{
        value=true //this is only activated to skip all the filtering - ie cellbender and restart with qc analysis once the parametes are changed
        file__anndata_merged = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/Deconvolution_Exp4_nfCore/results/merged_h5ad/adata.h5ad'
        file__cells_filtered = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/Deconvolution_Exp4_nfCore/results/merged_h5ad/adata-cell_filtered_per_experiment.tsv.gz'
    }
    
    run_celltype_assignment=true
    input_data_table = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/Deconvolution_Exp4_nfCore/inputs.tsv'
    run_with_genotype_input=true
	genotype_input {
        subset_genotypes = false
        full_vcf_file = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/nf_ci_deconv_inputs/donors_stephen_watts/lifted.vcf.gz'
    }
}
