params {

    run_gather_data = false // todo, check if used by pipeline
    
    input = 'cellbender' 
    extra_metadata = ''
    skip_preprocessing{
        value= false //this is only activated to skip all the filtering - ie cellbender and restart with qc analysis once the parametes are changed
        file__anndata_merged = ''
        file__cells_filtered = ''
    }
    
    run_celltype_assignment = true
    split_ad_per_bach=true //if not splitting the celltype assignment will be run on full tranche
    input_data_table = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/42603__ELGH_CardVal/inputs.tsv'
    run_with_genotype_input = false
	genotype_input {
        subset_genotypes = false
        full_vcf_file = ''
    }
}
