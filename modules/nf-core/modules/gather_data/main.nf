process GATHER_DATA{

    publishDir  path: "${outdir}",
                mode: "${params.copy_mode}",
                overwrite: "true"
    label 'process_medium'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/wtsihgi_nf_scrna_qc_0417190-2021-12-16-133460e8fb0b.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
    } else {
        container "wtsihgi/nf_scrna_qc:0417190"
    }

    input:
      path(outdir_prev)
      val(dummy_val)
      
    output:
      path("${subdir}", emit:outfiles_dataset)
      path("${subdir}_summary", emit:outfiles_dataset2)
      
      val(outdir, emit: outdir_dataset)

    script:
      outdir = "${outdir_prev}/handover"
      subdir = "minimal_dataset"
      if ("${params.input}" == 'cellranger'){
        cellbender_input='cellranger'
      }else if ("${params.input}" == 'existing_cellbender'){
        cellbender_input="${params.cellbender_file}"
      }else if("${params.input}" == 'cellbender'){
        cellbender_input='cellbender'
      }

      """
        gather_minimal_dataset.py \
          --output-dir=${subdir} \
          --results_dir=${outdir_prev} \
          --input_table=${params.input_data_table} \
          --cellbender=${cellbender_input} \
          --resolution=${params.resolution}
      """
}
