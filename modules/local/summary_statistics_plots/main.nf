process SUMMARY_STATISTICS_PLOTS {
        
    label 'process_low'
    publishDir  path: "${params.outdir}/handover",
            saveAs: { filename -> 
              filename == 'versions.yml' ? null : filename 
            },
            mode: "${params.copy_mode}",
            overwrite: "true"

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }
      
    input: 
        path(outdir_prev)
        val(gather_dummy_input)
        path(input_data_table)

    output: 
        path('Summary_plots'), emit: summary_plots
        path "versions.yml", emit: versions

    
    script:
      if ("${params.input}" == 'cellranger'){
          cellbender_input='cellranger'
        }else {
          cellbender_input='cellbender'
        }
        outdir = "${outdir_prev}/handover"
      """ 
          echo "${gather_dummy_input}" >dummy.out 
          transfer_data.py    --results_dir ${outdir_prev} \
                              --cb_res ${params.cellbender_resolution_to_use} \
                              --cellbender ${cellbender_input} \
                              --input_table ${input_data_table} \
                              --web_transfer ${params.webtransfer} \
                              --project_name ${params.project_name}
          cp ${params.extra_sample_metadata} Summary_plots/*/Summary || echo 'not available'
          cohort_report.py -d ${outdir_prev}
          
          cat <<-END_VERSIONS > versions.yml
          "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            python library argparse: \$(python -c "import argparse; print(argparse.__version__)")
            python library pandas: \$(python -c "import pandas; print(pandas.__version__)")
          END_VERSIONS

      """
}


process TRANSFER {
        
    label 'process_low'

    input: 
        path(summary_plots)
        path(rsync_to_web_file)
        path(results_dir)


    when:
        params.webtransfer
    script:

      """ 
        ./rsync_to_web.sh ${params.project_name} ${results_dir}       
      """
}