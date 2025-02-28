process DOUBLET_DECON{

    tag "${experiment_id}"
    // label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "mercury/azimuth_dsb:6_03_2024"
    }
    
    publishDir  path: "${params.outdir}/doublet_detection/DoubletDecon",
                mode: "${params.copy_mode}",
                overwrite: "true"
    memory { 
            sizeInGB = h5ad.size() / 1e9 * 5 * task.attempt
            return (sizeInGB ).toString() + 'GB' 
        }
        
    input:
        tuple(
            val(experiment_id),
            path(h5ad)
        )



    output:
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true
        tuple val(experiment_id), path("${experiment_id}__DoubletDecon_doublets_singlets.tsv"), emit: result  optional true
        // path("${experiment_id}__DoubletDetection_results.txt"), emit: doubletDetection_results

    script:
        
        outdir = "${params.outdir}/doublet_detection/"
        outdir = "${outdir}DoubletDecon"
        outfile = "${experiment_id}"

        """
            DoubletDecon.R -s ${h5ad} -o DoubletDecon_${experiment_id}
            
            if [ -e "DoubletDecon_${experiment_id}/DoubletDecon_doublets_singlets.tsv" ]; then
                ln -s DoubletDecon_${experiment_id}/DoubletDecon_doublets_singlets.tsv ${experiment_id}__DoubletDecon_doublets_singlets.tsv
            else
                echo "File DoubletDecon_${experiment_id}/DoubletDecon_doublets_singlets.tsv does not exist. Skipping symlink creation."
            fi
        """
}
