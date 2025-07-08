process SC_DBLFINDER {

    tag "${experiment_id}"
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }
    
    publishDir  path: "${params.outdir}/doublet_detection/scDblFinder",
                saveAs: { filename ->
                    if (filename.endsWith("scDblFinder_doublets_singlets.tsv")) {
                        return null
                    } else {
                        return filename
                    }
                },
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple(
            val(experiment_id),
            path(file_10x)
        )
    output:
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true
        tuple val(experiment_id), path("${experiment_id}__scDblFinder_doublets_singlets.tsv"), emit: result optional true
        path("scDblFinder_${experiment_id}")
    script:
        
        outdir = "${params.outdir}/"
        outdir = "${outdir}scDblFinder"
        outfile = "${experiment_id}"
        if (params.atac){
            atac = '--atac'
        }else{
            atac = ""
        }
        
        """
            export TMPDIR=\$PWD
            mkdir scDblFinder_${experiment_id}
            scDblFinder.R --tenX_matrix ${file_10x} -o scDblFinder_${experiment_id} ${atac}
            ln -s scDblFinder_${experiment_id}/scDblFinder_doublets_singlets.tsv ${experiment_id}__scDblFinder_doublets_singlets.tsv 
        """
}
