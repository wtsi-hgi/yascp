// Default parameters for harmony.
params.harmony = [
    run_process: false,
    variables_and_thetas: [value: [
        [variable: "experiment_id", theta: "1.0"],
        [variable: "experiment_id,phase", theta: "1.0,0.2"]
    ]]
]

def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

process HARMONY{
        // Takes PCs (rows = cell barcodes) and metadata (rows = cell barcodes),
    // runs Harmony
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    scratch false      // use tmp directory
    // storeDir '/tmp'
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("normalized_pca.h5ad")) {
                        null
                    } else if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else if(filename.endsWith("${param_details}.tsv.gz")) {
                        null
                    } else {
                        filename.replaceAll("${runid}-", "")
                    }
                },
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__metadata)
        path(file__pcs)
        val(params__pcs)
        each n_pcs
        each variables_and_thetas

    output:
        val(outdir, emit: outdir)
        path(file__anndata, emit: anndata)
        path(file__metadata, emit: metadata)
        path(file__pcs, emit: pcs)
        path("${runid}-reduced_dims.tsv.gz", emit: reduced_dims)
        // NOTE: passing the param details as an unpublished file is messy,
        // but I could not get collect of ${param_details} and file to work.
        path(
            "reduced_dims-${param_details}.tsv.gz",
            emit: reduced_dims_params
        )
        // val(n_pcs, emit: n_pcs)
        // tuple(
        //     val(outdir),
        //     path("reduced_dims.tsv.gz"),
        //     val(n_pcs),
        //     emit: results
        // )

    script:
        runid = "${params__pcs}_${n_pcs}_${variables_and_thetas.variable}"
        param_details = "${params__pcs}-harmony"
        param_details = "${param_details}.n_pcs=${n_pcs}"
        param_details = "${param_details}.variables=${variables_and_thetas.variable}"
        outdir = "${outdir_prev}/reduced_dims-${param_details}"
        outdir = "${outdir}.thetas=${variables_and_thetas.theta}"
        theta_str = "${variables_and_thetas.theta}".replaceAll("\\.", "pt")
        param_details = "${param_details}.thetas=${theta_str}"

        """
            0045-harmony_process_pcs.py \
                --pca_file ${file__pcs} \
                --metadata_file ${file__metadata} \
                --metadata_columns ${variables_and_thetas.variable} \
                --theta ${variables_and_thetas.theta} \
                --n_pcs ${n_pcs} \
                --out_file ${runid}-reduced_dims
            cp ${runid}-reduced_dims.tsv.gz \
            reduced_dims-${param_details}.tsv.gz
        """
        // NOTE: below code for harmony in R
        // 0045-harmony_process_pcs.R \
        //     --pca_file ${file__pcs} \
        //     --metadata_file ${file__metadata} \
        //     --metadata_columns ${variables_and_thetas.variable} \
        //     --theta ${variables_and_thetas.theta} \
        //     --n_pcs ${n_pcs} \
        //     --out_file reduced_dims
}
