process SC_DBLFINDER {

    tag "${experiment_id}"
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "mercury/azimuth_dsb:6_03_2024"
    }
    
    publishDir  path: "${params.outdir}/doublet_detection/scDblFinder",
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
        tuple val(experiment_id), path("${experiment_id}__scDblFinder_doublets_singlets.tsv"), emit: result optional true
        // path("${experiment_id}__DoubletDetection_results.txt"), emit: doubletDetection_results

    script:
        
        outdir = "${params.outdir}/"
        outdir = "${outdir}scDblFinder"
        outfile = "${experiment_id}"

        """

            mkdir TMP_DIR
            ln --physical ${file_10x_barcodes} TMP_DIR
            ln --physical ${file_10x_features} TMP_DIR
            ln --physical ${file_10x_matrix} TMP_DIR
            mkdir scDblFinder_${experiment_id}
            scDblFinder.R --tenX_matrix ./TMP_DIR --barcodes_filtered ${file_10x_barcodes} -o scDblFinder_${experiment_id}
            ln -s scDblFinder_${experiment_id}/scDblFinder_doublets_singlets.tsv ${experiment_id}__scDblFinder_doublets_singlets.tsv 
        """
}
