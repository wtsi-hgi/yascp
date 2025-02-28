process TOTAL_VI_INTEGRATION{
    
    if (params.utilise_gpu){
        label 'process_low'
    }else{
        label 'process_medium'
    }
    memory { 
            sizeInGB = adata.size() / 1e9 * 2 * task.attempt
            return (sizeInGB ).toString() + 'GB' 
        }

    publishDir  path: "${outdir_prev}/totalVi",
                mode: "${params.copy_mode}",
                overwrite: "true"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    input:
        path(adata)
        path(citedata)
        val(outdir_prev)

    output:
        path("./figures"), emit: figs, optional: true
        path("./scvi_model"), emit: scvi_model, optional: true
        path("./totalVI_integrated.h5ad"), emit: totalVI_integrated, optional: true
        
    script:

        """
            totalVI.py -h5ad_file ${adata} 
        """

}
