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
        path('Summary_plots')
    
    script:
      if ("${params.input}" == 'cellranger'){
          cellbender_input='cellranger'
        }else if ("${params.input}" == 'existing_cellbender'){
          cellbender_input="${params.cellbender_file}"
        }else if("${params.input}" == 'cellbender'){
          cellbender_input='cellbender'
        }
        outdir = "${outdir_prev}/handover"
      """
          transfer_data.py    --results_dir ${outdir_prev} \
                              --cb_res ${params.resolution} \
                              --cellbender ${cellbender_input} \
                              --input_table ${input_data_table}
          
      """
}
