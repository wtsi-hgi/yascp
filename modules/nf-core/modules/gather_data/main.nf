process GATHER_DATA{

    publishDir  path: "${outdir}",
                mode: "${params.copy_mode}",
                overwrite: "true"
    label 'process_medium'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }

    input:
      path(outdir_prev)
      
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