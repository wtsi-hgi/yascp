process metadata_posthoc
{

  label 'process_tiny'

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "/software/hgi/containers/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
      //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
  } else {
      container "wtsihgi/nf_scrna_qc:6bb6af5"
  }

  input:
    tuple(val(exp_id),path(output_dir))

  output:
    path("out.txt"), emit: dummy_out

  script:

    """
      update_metadata_posthoc.py --exp_id ${exp_id} --res ${output_dir} --met_folder ${params.extra_metadata}
      cp results/handover/Donor_Quantification_summary/${exp_id}_Donor_Report.tsv results/handover/Summary_plots/${exp_id}/Summary/${exp_id}_Donor_Report.tsv
      cp results/handover/Donor_Quantification_summary/${exp_id}_Tranche_Report.tsv results/handover/Summary_plots/${exp_id}/Summary/${exp_id}_Tranche_Report.tsv
      echo Done > out.txt
    """
}

process replace_donors_posthoc
{


  label 'process_tiny'

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      container "/software/hgi/containers/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
      //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
  } else {
      container "wtsihgi/nf_scrna_qc:6bb6af5"
  }

  input:
    tuple(val(exp_id),path(results))

  output:
    path("out.txt"), emit: dummy_out

  script:

    """

      replace_donors_posthoc.py -i ${results}/deconvolution/vireo_gt_fix/assignments_all_pools.tsv --genotype_phenotype_mapping ${params.genotype_phenotype_mapping_file} --input_file "results/handover/Summary_plots/${exp_id}/Fetch Pipeline/Input/input_table.tsv"
      echo Done > out.txt
    """
}