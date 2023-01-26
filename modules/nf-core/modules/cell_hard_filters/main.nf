process CELL_HARD_FILTERS{
    tag "${samplename}"
    
    label 'process_medium_single_CPU' 
    label 'process_medium_memory'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    input:
        path(file_paths_h5ad)
        val(hard_filter_file_params)
        val(drop)

    // NOTE: use path here and not file see:
    //       https://github.com/nextflow-io/nextflow/issues/1414
    when:
        hard_filter_file_params != "no_file__file_sample_qc"
    output:
        path("hard_filters_${file_paths_h5ad}", emit: anndata)

    script:
        """
        cell_hard_filters.py --params_yaml ${hard_filter_file_params} \
            --metadata_key ${params.metadata_key_column.value} \
            --number_cpu ${task.cpus} \
            --output_file hard_filters_adata \
            --h5addata_file ${file_paths_h5ad} \
            --anndata_compression_opts ${params.anndata_compression_opts} \
            --drop ${drop}
        """

}