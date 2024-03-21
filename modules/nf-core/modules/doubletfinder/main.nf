process DOUBLET_FINDER {

    tag "${experiment_id}"
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/azimuth_dsb_6_03_2024.sif"
    } else {
        container "mercury/azimuth_dsb:6_03_2024"
    }
    
    publishDir  path: "${params.outdir}/multiplet.method=DoubletFinder",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple(
            val(experiment_id),
            path(h5ad)
        )
        val(expected_multiplet_rate)

    output:
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true
        tuple val(experiment_id), path("${experiment_id}__DoubletFinder_doublets_singlets.tsv"), emit: result
        // path("${experiment_id}__DoubletDetection_results.txt"), emit: doubletDetection_results

    script:
        
        outdir = "${params.outdir}/multiplet"
        outdir = "${outdir}.method=doubletFinder"
        outfile = "${experiment_id}"

        """
            DoubletFinder.R -s ${h5ad} -c FALSE -o DoubletFinder_${experiment_id} --expected_multiplet_rate ${expected_multiplet_rate}
            ln -s DoubletFinder_${experiment_id}/DoubletFinder_doublets_singlets.tsv ${experiment_id}__DoubletFinder_doublets_singlets.tsv
        """
}
