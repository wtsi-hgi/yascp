process TOTAL_VI_INTEGRATION{
    
    if (params.utilise_gpu){
        label 'process_low'
    }else{
        label 'process_low'
    }

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/sle_project/5.sle_noCB_fullRun_MAD_perPool/yascp_totalvi_v1.sif"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    input:
        path(adata)
        path(citedata)

    // output:
    //     val(outdir, emit: outdir)
        
    script:

        """
            totalVI.py -h5ad_file ${adata} 
        """

}
