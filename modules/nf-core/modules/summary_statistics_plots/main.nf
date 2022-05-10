process SUMMARY_STATISTICS_PLOTS {
        
    label 'process_low'
    publishDir  path: "${outdir}",
            mode: "${params.copy_mode}",
            overwrite: "true"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }

    input: 
        path(outdir_prev)
        val(gather_dummy_input)
        path(input_data_table)

    output: 
        path('Summary_plots'), emit: summary_plots
    
    script:
      if ("${params.input}" == 'cellranger'){
          cellbender_input='cellranger'
        }else {
          cellbender_input='cellbender'
        }
        outdir = "${outdir_prev}/handover"
      """
          transfer_data.py    --results_dir ${outdir_prev} \
                              --cb_res ${params.cellbender_resolution_to_use} \
                              --cellbender ${cellbender_input} \
                              --input_table ${input_data_table} \
                              --web_transfer ${params.webtransfer} \
                              --project_name ${params.project_name}
          
      """
}


process TRANSFER {
        
    label 'process_low'

    input: 
        path(summary_plots)


    when:
        params.webtransfer
    script:

      """
        ../../../scripts/rsync_to_web.sh ${params.project_name}          
      """
}