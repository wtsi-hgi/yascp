
params {
    input = 'cellbender'
    cellbender_file='' //if cellbender is run already then can skip this by selecting existing_cellbender and input 
    extra_metadata = ''
    skip_preprocessing{
        value=false //this is only activated to skip all the filtering - ie cellbender and restart with qc analysis once the parametes are changed
        file__anndata_merged = ''
        file__cells_filtered = ''
    }
    
    run_celltype_assignment=true
    input_data_table = '../inputs.tsv'
    run_with_genotype_input=true
	genotype_input {
        subset_genotypes = false
        full_vcf_file = 'lifted.vcf.gz'
    }
}

