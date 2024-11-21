process SCDS {

    tag "${experiment_id}"
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.azimuth_dsb_container}"
    } else {
        container "mercury/azimuth_dsb:6_03_2024"
    }
    
    publishDir  path: "${params.outdir}/doublet_detection/multiplet.method=SCDS",
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
        tuple val(experiment_id), path("${experiment_id}__scds_doublets_singlets.tsv"), emit: result
        // path("${experiment_id}__DoubletDetection_results.txt"), emit: doubletDetection_results

    script:
        
        outdir = "${params.outdir}/multiplet"
        outdir = "${outdir}.method=SCDS"
        outfile = "${experiment_id}"

        """

            mkdir TMP_DIR
            ln --physical ${file_10x_barcodes} TMP_DIR
            ln --physical ${file_10x_features} TMP_DIR
            ln --physical ${file_10x_matrix} TMP_DIR
            mkdir scds_${experiment_id}
            scds.R -t ./TMP_DIR -o scds_${experiment_id}
            cat scds_${experiment_id}/scds_doublets_singlets.tsv | awk -F' ' '{print \$1"\\t"\$3"\\t"\$2}' > ${experiment_id}__scds_doublets_singlets.tsv
        """
}