process TOTAL_VI_INTEGRATION{
    
    if (params.utilise_gpu){
        label 'process_low'
    }else{
        label 'process_medium'
    }

    publishDir  path: "${outdir_prev}/clustering/totalVi",
                mode: "${params.copy_mode}",
                overwrite: "true"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/yascp_totalvi_v1.sif"
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
