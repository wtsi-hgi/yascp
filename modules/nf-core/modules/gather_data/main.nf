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

    when:
      params.data_handover.run_process

    input:
      val(outdir_prev)
      val(cellbender_input)

    output:
      path("${subdir}/*", emit:outfiles_dataset)
      path("${subdir2}/*", emit:outfiles_dataset2)
      val(outdir, emit: outdir_dataset)

    script:
      outdir = "${outdir_prev}/handover"
      subdir = "minimal_dataset"
      subdir2="${subdir}_summary"
      if (multiplet_calls) {
        argstr = " --scrublet-output-dir=${qc_output_dir}/multiplet.method=scrublet"
      } else {
        argstr = ""
      }
      """
        gather_minimal_dataset.py \
          --output-dir=${subdir} \
          --results_dir= ${outdir_prev} \
          --input_table= ${params.input_data_table} \
          --cellbender = ${cellbender_input}
      """
}