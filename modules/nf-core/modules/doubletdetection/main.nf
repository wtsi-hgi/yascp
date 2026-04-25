process DOUBLET_DETECTION {

    tag "${experiment_id}"
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.nf_scrna_qc_v3_container}"
    } else {
        container "mercury/nf_scrna_qc:v3"
    }
    
    publishDir  path: "${params.outdir}/doublet_detection/multiplet.method=doubletdetection",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple(
            val(experiment_id),
            path(file_10x_barcodes),
            path(file_10x_features),
            path(file_10x_matrix)
        )

    output:
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true
        tuple val(experiment_id), path("${experiment_id}__DoubletDetection_results.txt"), emit: result

    script:
        
        outdir = "${params.outdir}/doublet_detection/multiplet"
        outdir = "${outdir}.method=doubletdetection"
        outfile = "${experiment_id}"

        """
            mkdir TMP_DIR
            ln --physical ${file_10x_barcodes} TMP_DIR
            ln --physical ${file_10x_features} TMP_DIR
            ln --physical ${file_10x_matrix} TMP_DIR
            DoubletDetection.py --tenxdata_dir ./TMP_DIR --n_iterations 100
            ln -s DoubletDetection_results.txt ${experiment_id}__DoubletDetection_results.txt
        """
}