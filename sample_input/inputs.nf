
params {
    extra_metadata = '/path/to/extra_metadata.tsv'   //Sometimes users may want to merge extra known metadata for a pool in the h5ad files prior to qc
    extra_sample_metadata ="/path/to/donor_extra_metadata.tsv"  //Sometimes users may want to merge extra known metadata for a donor within pool prior to qc
    input_data_table = '/full/path/to/input_table.tsv' //Required!! This points to all the cellranger files and pool definition files.
    genotype_input {
        run_with_genotype_input=true //if false do not need the genotype_input parameters.
        vireo_with_gt=false // Vireo is capable in runing both with genotypes and without. Here we define in which mode we want to run it.
        posterior_assignment = false //if this is set to true, we will perform the genotype donor matching after the deconvolution is performed.
        subset_genotypes = false
        tsv_donor_panel_vcfs = "/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/OneK1k/vcf_inputs.tsv" //this is a panel of vcf files that we want to compar the genotypes with
    }
}

