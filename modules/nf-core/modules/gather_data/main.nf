process GATHER_DATA{

    publishDir  path: "${params.outdir}/handover",
                mode: "${params.copy_mode}",
                overwrite: "true"
    label 'process_medium'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
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
      path("Donor_Quantification/*/*.tsv", emit: barcodes_files)
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
        echo ${dummy_val}
        gather_minimal_dataset.py \
          --output-dir=${subdir} \
          --results_dir=${outdir_prev} \
          --input_table=${input_data_table} \
          --cellbender=${cellbender_input} \
          --resolution=${params.cellbender_resolution_to_use} \
          --write_h5=${params.write_h5} \
          --experiment_name=${params.RUN}
      """
}

process SPLIT_DATA_BY_STUDY
{
  publishDir  path: "${params.outdir}",
              mode: "${params.copy_mode}",
              overwrite: "true"

  label 'process_tiny'

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
      //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
  } else {
      container "wtsihgi/nf_scrna_qc:6bb6af5"
  }

  input:
    path(outdir_prev)
    path(input_dir)
    path(input_dirs_cram)
    val(dummy_val)

  output:
    path("${outdir_ukbb}", emit: outdir_ukbb)

  script:
    donor_assignments_tsv = "${outdir_prev}/deconvolution/vireo_gt_fix/assignments_all_pools.tsv"
    outdir = "${outdir_prev}"
    outdir_ukbb = "handover_study/GT_UKBB"
    """
      echo "${input_dirs_cram}" > ./cram_dirs.lst
      split_dataset_by_study.py ${donor_assignments_tsv} ${input_dir} ./handover_study ./cram_dirs.lst
    """
}
