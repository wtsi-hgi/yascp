process VALIDATE_YASCP_OUTPUTS{
    publishDir  "${params.outdir}/validation_test/", mode: 'copy'
    label 'process_tiny'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

  input:
    //tuple val(pool), path (input_file)
    path input_file
    val azimuth_run
    val celltypist_run
    val scpred_run
    //tuple val(step), path (celltype_validation)
    //tuple val(pool), path (vireo_validation)
    //tuple val(pool), path (gt_fix_header_validation)
    //tuple val(pool), path (cellsnp_validation)
    //path(genotype_matcher_validation)
    //tuple val(pool), path (preprocessing_citeseq_validation)
    //tuple val(pool), path (preprocessing_citeseq_filtered_validation)
    //path(pca_validation)

  output:
      path "validation_data_combined.tsv", emit: validation_res
      path "results_validation_report.txt", emit: validation_rep
    script:
  """
    validate_output.py ${azimuth_run} ${celltypist_run} ${scpred_run}

  """
}

process RENAME_VALIDATION_FILES {
    input:
    tuple val(sample_id), path(input_file)
    
    output:
    tuple val(sample_id), path("*_filtered.counts.txt"), emit:output_check
    
   
    script:
    def prefix = input_file.name.tokenize('.')[0]
    """
    #cp ${input_file} ${prefix}_filtered.counts.txt
    ln -s \$(readlink -f ${input_file}) ${prefix}_filtered.counts.txt
    """
}

process GET_CELLBENDER_NUMBERS {
    input:
    tuple val(sample_id), path(input_folder)
    
    output:
    tuple val(sample_id), path("${sample_id}_cellbender.counts.txt"), emit:output_check
    
   
    script:
    """
    if [[ -f ${input_folder}/barcodes.tsv.gz ]]; then
      zcat ${input_folder}/barcodes.tsv.gz
    else
      cat ${input_folder}/barcodes.tsv
    fi | wc -l > ${sample_id}_cellbender.counts.txt
    """
}