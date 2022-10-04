
process GT_IBD_MATCH_RELATIONSHIPS
{
  tag "${pool_panel_id}"
  publishDir  path: "${params.outdir}/gtmatch/${pool_id}",
          pattern: "*.csv",
          mode: "${params.copy_mode}",
          overwrite: "true"

  if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
      // println "container: /software/hgi/containers/wtsihgi-nf_genotype_match-1.0.sif\n"
      container "/software/hgi/containers/wtsihgi-nf_genotype_match-1.0.sif"
  } else {
      container "mercury/wtsihgi-nf_genotype_match-1.0"
  }

  input:
    tuple val(pool_panel_id), path(gtcheck_output_files)

  output:
    
  label 'process_medium'

  script:

  """
    echo 'lets do this - performing the relationship check with yielded matches'
    
  """
}