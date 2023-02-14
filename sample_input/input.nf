
params {
    cellbender_file='' //if cellbender is run already then can skip this by selecting existing_cellbender and input 
    extra_metadata = '/path/to/extra_metadata.tsv'   
    input_data_table = '/full/path/to/input_table.tsv'
    
	genotype_input {
        run_with_genotype_input=true //if false do not need the genotype_input parameters.
        subset_genotypes = false
        full_vcf_file = 'lifted.vcf.gz' //vcf file in hg38 format without chr prefix
        posterior_assignment = false //this allows running vireo and assign donors posthoc.
    }
}

