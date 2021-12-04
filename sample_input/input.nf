params {

    extra_metadata = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/mo11_work/nfCore/standalone/mo11_scRNA_dev/sample_input/extra_metadata.txt'
    skip_preprocessing{
        value=false //this is only activated to skip all the filtering - ie cellbender and restart with qc analysis once the parametes are changed
        file__anndata_merged = '/lustre/scratch119/humgen/projects/anderson_organoids/otar_mucosal/stimulation_IFNg/full_run/results/merged_h5ad/adata.h5ad'
        file__cells_filtered = '/lustre/scratch119/humgen/projects/anderson_organoids/otar_mucosal/stimulation_IFNg/full_run/results/merged_h5ad/adata-cell_filtered_per_experiment.tsv.gz'
    }
    input = 'cellranger'
    input_data_table = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/mo11_work/nfCore/standalone/mo11_scRNA_dev/sample_input/input_table.tsv'
    run_with_genotype_input=false
	genotype_input {
        subset_genotypes = false
        full_vcf_file = ''
    }

}
