process GATHER_DATA{

    publishDir  path: "${outdir}",
                mode: "${params.copy_mode}",
                overwrite: "true"
    label 'process_medium'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    input:
      path(outdir_prev)
      val(dummy_val)
      path(input_data_table)
    output:
      path("${subdir}", emit:outfiles_dataset)
      path("${subdir}_summary", emit:outfiles_dataset2)

      val(outdir, emit: outdir_dataset)

    script:
      outdir = "${outdir_prev}/handover"
      subdir = "Donor_Quantification"
      if ("${params.input}" == 'cellranger'){
        cellbender_input='cellranger'
      }else{
        cellbender_input='cellbender'
      }

      """
        gather_minimal_dataset.py \
          --output-dir=${subdir} \
          --results_dir=${outdir_prev} \
          --input_table=${input_data_table} \
          --cellbender=${cellbender_input} \
          --resolution=${params.cellbender_resolution_to_use} \
          --write_h5=${params.write_h5} \

      """
}

process SPLIT_DATA_BY_STUDY
{
  publishDir  path: "${outdir}",
              mode: "${params.copy_mode}",
              overwrite: "true"

  label 'process_small'

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "/software/hgi/containers/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
      //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
  } else {
      container "wtsihgi/nf_scrna_qc:6bb6af5"
  }

  input:
    path(outdir_prev)
    path(input_dir)
    tuple val(pool_id), path(donor_assignments_csv)

  output:
    path("${outdir_ukbb}", emit: outdir_ukbb)

  script:
    outdir = "${outdir_prev}/handover_study"
    outdir_ukbb = "${outdir}/GT_UKBB"
    """
      split_dataset_by_study.py ${donor_assignments_csv} ${input_dir} ${pool_id}
    """
}
