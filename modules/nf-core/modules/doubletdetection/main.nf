process DOUBLET_DETECTION {

    tag "${experiment_id}"
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/nf_scrna_qc_v3.img"
    } else {
        container "mercury/nf_scrna_qc:v3"
    }
    
    publishDir  path: "${params.outdir}/multiplet.method=doubletdetection",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        // path(outdir_prev)
        tuple(
            val(experiment_id),
            path(file_10x_barcodes),
            path(file_10x_features),
            path(file_10x_matrix)
        )

    // output:
    //     path("plots/*.pdf") optional true
    //     path("plots/*.png") optional true

    script:
        
        outdir = "${params.outdir}/multiplet"
        outdir = "${outdir}.method=doubletdetection"
        outfile = "${experiment_id}"

        """
            TMP_DIR=\$(mktemp -d -p \$(pwd))
            ln --physical ${file_10x_barcodes} \$TMP_DIR
            ln --physical ${file_10x_features} \$TMP_DIR
            ln --physical ${file_10x_matrix} \$TMP_DIR
            DoubletDetection.py --tenxdata_dir \$TMP_DIR --n_iterations 100
        """
}