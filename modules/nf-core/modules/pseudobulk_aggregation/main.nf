process PSEUDOBULK_AGREGATION {
  publishDir  path: "${outdir}/aggregated_counts",mode: "${params.copy_mode}",
              overwrite: "true"
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
        
    } else {
        container "${params.yascp_container_docker}"
    }

  input:
    path(adata) // lists input files per donor
    val(agg_column)
    val(n_cells_min)
    val(n_donors_min)
    path(assignments_all_pools)

  output:
    path("*phenotype_file.tsv", emit:phenotype_file)
    path("*genotype_phenotype_mapping.tsv", emit:genotype_phenotype_mapping)
    path('*.tsv')

  script:
  outdir = params.outdir
  if ("${assignments_all_pools}" == 'assignments_all_pools.tsv'){
    generate_gp_mapping = 'preprocess_gt_assignments.py --infile assignments_all_pools.tsv'
    agg_couns = "--genotype_phenotype genotype_phenotype.tsv"
  }else{
    generate_gp_mapping = ''
    agg_couns = ""
  }
  """
    ${generate_gp_mapping}
    echo ${assignments_all_pools}
    echo 'Lets aggregate the cells in the andata: ${adata} object is pseudobulk'
    aggregate_sc_data.py --agg_column '${agg_column}' --n_cells ${n_cells_min} -n_individ ${n_donors_min} -h5ad ${adata} --method ${params.eQTL.aggregation_method} ${agg_couns}
  """
}