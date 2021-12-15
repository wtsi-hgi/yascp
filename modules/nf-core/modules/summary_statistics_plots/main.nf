process SUMMARY_STATISTICS_PLOTS {
        
    label 'process_low'
    publishDir  path: "${outdir}",
            mode: "${params.copy_mode}",
            overwrite: "true"
    
    

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/scrna_deconvolution_latest.img"
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
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
