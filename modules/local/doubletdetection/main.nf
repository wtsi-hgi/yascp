process DOUBLET_DETECTION {

    tag "${experiment_id}"
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }
    
    publishDir  path: "${params.outdir}/doublet_detection/multiplet.method=doubletdetection",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple(
            val(experiment_id),
            path(gex_h5ad)
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
            DoubletDetection.py --tenxdata_dir ${gex_h5ad} --n_iterations 100
            ln -s DoubletDetection_results.txt ${experiment_id}__DoubletDetection_results.txt
        """
}