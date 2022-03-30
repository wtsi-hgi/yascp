
params {
    input = 'cellbender'
    cellbender_file='' //if cellbender is run already then can skip this by selecting existing_cellbender and input 
    extra_metadata = '/path/to/extra_metadate.tsv'   
    input_data_table = '/full/path/to/inputs.tsv'
    run_with_genotype_input=true //if false do not need the genotype_input parameters.
	genotype_input {
        subset_genotypes = false
        full_vcf_file = 'lifted.vcf.gz'
        posterior_assignment = false //this allows running vireo and assign donors pothoc.
    }
}

