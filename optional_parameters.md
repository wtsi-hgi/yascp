# Optional Parameters
Multiple required/optional inputs are described below. 

```console
params {
    //REQUIRED
    input_data_table = '/path/to/input.tsv' //A samplesheet file containing paths to all the cellranger and pool definition files

    //OPTIONAL
    split_ad_per_bach=true //This parameter defines whether cell type assignment is run on the full dataset together (false) or per batch (true)

    extra_metadata = '/path/to/extra_metadata.tsv' //A file with extra known metadata to merge for a pool in the h5ad files prior to QC

    extra_sample_metadata ='/path/to/donor_extra_metadata.tsv' //A file with extra known metadata to merge for a donor within a pool prior to QC

    //cellbender_location='/path/to/existing/folder/nf-preprocessing/cellbender' //!!!!! Uncomment this and edit the path, if cellbender results are already available (even partial results). The pipeline will skip the cellbender step for samples that already have cellbender results.

    existing_cellsnp="" //Provide a path to cellsnp results (if they are already available, even partial results) to skip cellsnp step for the files with results.

    genotype_input {
        run_with_genotype_input=true //This parameter defines whether the genotype_input is used (true) or not(false). If this is set to true tsv_donor_panel_vcfs has to be specified
        tsv_donor_panel_vcfs = "/path/to/reference/panel/vcf_inputs.tsv" //A file containing paths to vcf files with a priori known genotypes that we want to compare the genotypes from samples with
        vireo_with_gt=false //This parameter defines whether Vireo is run with a priori known genotypes (true) or not (false)
        posterior_assignment = false //if this is set to true, and a priori known genotypes are provided, after deconvolution the genotypes will be matched to Vireo-detected donors
        subset_genotypes = false //This parameter defines whether to subset a large genotype file to include only the genotypes expected in a pool to reduce the deconvolution time (true) or not (false).
    }
}


```
### Required parameters
`input_data_table` - a samplesheet file containing paths to all the cellranger and pool definition files.

### Optional parameters
`split_ad_per_bach` - this parameter defines whether cell type assignment is run on the full dataset together (false) or per batch (true).

`extra_metadata` - a file with extra known metadata to merge for a pool in the h5ad files prior to QC.

`extra_sample_metadata` - a file with extra known metadata to merge for a donor within a pool prior to QC.

`cellbender_location` - uncomment this and edit the path, if cellbender results are already available (even partial results). The pipeline will skip the cellbender step for samples that already have cellbender results. For more details see `Some tricks to avoid rerunning the pipeline over and over if you already have some partial data`

`existing_cellsnp` - provide a path to cellsnp results (if they are already available, even partial results) to skip cellsnp step for the files whth results. For more details see `Some tricks to avoid rerunning the pipeline over and over if you already have some partial data`

`run_with_genotype_input` - this parameter defines whether the genotype_input is used (true) or not(false). If this is set to true tsv_donor_panel_vcfs has to be specified.

`tsv_donor_panel_vcfs` - a file containing paths to vcf files with a priori known genotypes that we want to compare the genotypes from samples with.

`vireo_with_gt` - this parameter defines whether Vireo is run with a priori known genotypes (true) or not (false).

`posterior_assignment` if this is set to true, and a priori known genotypes are provided, after deconvolution the genotypes will be matched to Vireo-detected donors

`subset_genotypes` - this parameter defines whether to subset a large genotype file to include only the genotypes expected in a pool to reduce the deconvolution time (true) or not (false).
