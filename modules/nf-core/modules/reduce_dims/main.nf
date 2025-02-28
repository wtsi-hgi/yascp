
process DONT_INTEGRATE{

    scratch false      // use tmp directory
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"       
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__metadata)
        path(file__pcs)
        val(params__pcs)
        each n_pcs

    output:
        val(outdir, emit: outdir)
        path("reduced_dims.tsv.gz", emit: reduced_dims)
        path("${outfile}-dontIntegrate.h5ad", emit: anndata)
        // NOTE: passing the param details as an unpublished file is messy,
        // but I could not get collect of ${param_details} and file to work.
        path(
            "reduced_dims-${param_details}.tsv.gz",
            emit: reduced_dims_params
        )

    script:
        param_details = "${params__pcs}-dontIntegrate"
        param_details = "${param_details}.n_pcs=${n_pcs}"
        outdir = "${outdir_prev}/reduced_dims-${param_details}"
        outfile = "outfile_adata"
        """
        
        0045-reduce_dims_file.py \
            --h5_anndata ${file__anndata} \
            --n_pcs ${n_pcs} \
            --output_file ${outfile}-dontIntegrate
        cp ${outfile}-dontIntegrate-reduced_dims.tsv.gz \
            reduced_dims-${param_details}.tsv.gz
        mv ${outfile}-dontIntegrate-reduced_dims.tsv.gz \
            reduced_dims.tsv.gz
        """

}
