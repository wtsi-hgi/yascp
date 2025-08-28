process CELL_HARD_FILTERS{
    tag "${samplename}"
    label 'process_high'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }


    publishDir  path: "${params.outdir}/handover/merged_h5ad/",
                saveAs: {filename ->
                    if (filename.contains("hard_filters_")) {
                        filename = '2.hard_filters_annotated_h5ad.h5ad'
                    }else{
                        filename
                    }
                },
                mode: "${params.copy_mode}",
                overwrite: "true"


    input:
        path(file_paths_h5ad)
        val(drop)

    // NOTE: use path here and not file see:
    //       https://github.com/nextflow-io/nextflow/issues/1414
    when:
        params.sample_qc.cell_filters.experiment.value != '' | params.sample_qc.cell_filters.all_samples.value != '' | params.sample_qc.downsample_cells_fraction.value != '' | params.sample_qc.downsample_cells_n.value != '' | params.sample_qc.downsample_feature_counts.value != ''
    output:
        path("hard_filters_*.h5ad", emit: anndata)

    script:


        if (params.sample_qc.cell_filters.experiment.value != ''){
            cell_filters_experiment =  "--cell_filters_experiment \"${params.sample_qc.cell_filters.experiment.value}\""
        }else{
            cell_filters_experiment =  ""
        }

        if (params.sample_qc.cell_filters.all_samples.value != ''){
            cell_filters =  "--cell_filters '${params.sample_qc.cell_filters.all_samples.value}'"
        }else{
            cell_filters =  ""
        }

        if (params.sample_qc.downsample_cells_fraction.value != ''){
            downsample_cells_fraction = "--downsample_cells_fraction ${params.sample_qc.downsample_cells_fraction.value}"
        }else{
            downsample_cells_fraction =  ""
        }


        if (params.sample_qc.downsample_cells_n.value != ''){
            downsample_cells_n =  "--downsample_cells_n ${params.sample_qc.downsample_cells_n.value}"
        }else{
            downsample_cells_n =  ""
        }

        if (params.sample_qc.downsample_feature_counts.value != ''){
            downsample_feature_counts =  "--downsample_feature_counts ${params.sample_qc.downsample_feature_counts.value}"
        }else{
            downsample_feature_counts =  ""
        }
        
        """
        cell_hard_filters.py \
            ${downsample_cells_fraction} \
            ${cell_filters} \
            ${downsample_cells_n} \
            ${downsample_feature_counts} \
            ${cell_filters_experiment} \
            --metadata_key ${params.metadata_key_column.value} \
            --number_cpu ${task.cpus} \
            --output_file hard_filters_adata \
            --h5addata_file ${file_paths_h5ad} \
            --anndata_compression_opts ${params.anndata_compression_opts} \
            --drop ${drop}
        """

}