def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

process dummy_filtered_channel{
    label 'process_low' 
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    input:
        path(file_paths_h5ad)
        val(id_in)
    output:
        path("filtered_cell_dummy.tsv", emit: anndata_metadata)
        
    script:
        """
            echo 'lets do it' > test.txt
            dummy_filtered_channel.py --h5_anndata ${file_paths_h5ad} -id ${id_in}
        """

}




process merge_samples_from_h5ad {
    // Takes a list of h5ad files and merges them into one anndata object.
    // ------------------------------------------------------------------------
    //cache true        // cache results from run
    tag "${samplename}"
    
    label 'process_medium_single_CPU' 
    label 'process_medium_memory'

    publishDir  path: "${outdir}/handover/merged_h5ad",
                saveAs: {filename -> filename.replaceAll("-", "pre_QC_")},
                mode: "${params.copy_mode}",
                overwrite: "true"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    input:
        val(outdir_prev)
        path(file_paths_h5ad)
        path(file_metadata)
        val(file_cellmetadata)
        val(metadata_key)
        file(file_h5ad)
        val(anndata_compression_opts)
        path(celltype)

    // NOTE: use path here and not file see:
    //       https://github.com/nextflow-io/nextflow/issues/1414
    output:
        path("1.pre_QC_adata.h5ad", emit: anndata)
        // path(
        //     "pre_QC_adata-cell_filtered_per_experiment.tsv.gz",
        //     emit: cells_filtered
        // )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        outdir = "${outdir_prev}"
        // String filename = './parameters.yml'
        // yaml.dump(file_params , new FileWriter(filename))
        // Customize command for optional files.
        if (params.extra_metadata!=''){
            extra_metadata = "--extra_metadata ${params.extra_metadata}"
        }else{
            extra_metadata = ""
        }

        if (params.extra_sample_metadata!='' && params.add_donor_metadata){
            // If same columns defined in the h5ad, these will be overwritten
            // If same columns in extra donor metadata as in extra sample metadata, sample metadata will be overwritten
            extra_sample_metadata = "add_extra_sample_metadata.py --vireo ${file_metadata} --extra_sample_metadata ${params.extra_sample_metadata} --metadata_key ${metadata_key}"
        }else{
            extra_sample_metadata = "ln -s ${file_metadata} Vireo_metadata.csv"
        }


        cmd__params = ""
        cmd__cellmetadata = ""
        if (file_cellmetadata != "no_file__file_cellmetadata") {
            cmd__cellmetadata = "--cell_metadata_file ${file_cellmetadata}"
        }
        files__h5ad = file_h5ad.join(',')
        """
        echo "publish_directory: ${celltype}"
        rm -fr plots
        nf_helper__prep_h5addata_file.py \
            --h5ad_list ${files__h5ad} \
            --output_file nf_prepped__file_paths_h5ad.tsv
        ${extra_sample_metadata}
        scanpy_merge_from_h5ad.py \
            --h5addata_file nf_prepped__file_paths_h5ad.tsv \
            --sample_metadata_file Vireo_metadata.csv \
            --sample_metadata_columns_delete "sample_status,study,study_id" \
            --metadata_key ${metadata_key} \
            --number_cpu ${task.cpus} \
            --output_file 1.pre_QC_adata \
            --anndata_compression_opts ${anndata_compression_opts} --celltype ${celltype} \
            ${cmd__params} \
            ${cmd__cellmetadata} ${extra_metadata}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}

process merge_samples {
    // Takes a list of raw 10x files and merges them into one anndata object.
    // ------------------------------------------------------------------------
    //cache true        // cache results from run
    
    tag "${samplename}"
    
    label 'process_high'
    publishDir  path: "${outdir}/handover/merged_h5ad",
                saveAs: {filename -> filename.replaceAll("-", "pre_QC_")},
                mode: "${params.copy_mode}",
                overwrite: "true"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    input:
        val(outdir_prev)
        path(file_paths_10x)
        path(file_metadata)
        val(file_cellmetadata)
        val(metadata_key)
        file(file_10x_barcodes)
        file(file_10x_features)
        file(file_10x_matrix)
        val(anndata_compression_opts)

    // NOTE: use path here and not file see:
    //       https://github.com/nextflow-io/nextflow/issues/1414
    output:
        path("adata.h5ad", emit: anndata)
        // path(
        //     "adata-cell_filtered_per_experiment.tsv.gz",
        //     emit: cells_filtered
        // )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:

        outdir = "${outdir_prev}"
        // String filename = './parameters.yml'
        // yaml.dump(file_params , new FileWriter(filename))
        // Customize command for optional files.
        cmd__params = ""
        cmd__cellmetadata = ""
        if (file_cellmetadata != "no_file__file_cellmetadata") {
            cmd__cellmetadata = "--cell_metadata_file ${file_cellmetadata}"
        }
        files__barcodes = file_10x_barcodes.join(',')
        files__features = file_10x_features.join(',')
        files__matrix = file_10x_matrix.join(',')

        """
        echo "publish_directory: ${outdir}"
        rm -fr plots
        nf_helper__prep_tenxdata_file.py \
            --barcodes_list ${files__barcodes} \
            --features_list ${files__features} \
            --matrix_list ${files__matrix} \
            --tenxdata_file ${file_paths_10x} \
            --output_file nf_prepped__file_paths_10x.tsv
        scanpy_merge.py \
            --tenxdata_file nf_prepped__file_paths_10x.tsv \
            --sample_metadata_file ${file_metadata} \
            --sample_metadata_columns_delete "sample_status,study,study_id" \
            --metadata_key ${metadata_key} \
            --number_cpu ${task.cpus} \
            --output_file adata \
            --anndata_compression_opts ${anndata_compression_opts} \
            ${cmd__params} \
            ${cmd__cellmetadata}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process prep_merge_samples {
    input:
        tuple(
            val(experiment_id),
            path(file_10x_barcodes),
            path(file_10x_features),
            path(file_10x_matrix)
        )
    label 'process_tiny'
    output:
        path("${experiment_id}---barcodes.tsv.gz", emit: barcodes)
        path("${experiment_id}---features.tsv.gz", emit: features)
        path("${experiment_id}---matrix.mtx.gz", emit: matrix)

    script:
        """
        ln --physical ${file_10x_barcodes} ${experiment_id}---barcodes.tsv.gz
        ln --physical ${file_10x_features} ${experiment_id}---features.tsv.gz
        ln --physical ${file_10x_matrix} ${experiment_id}---matrix.mtx.gz
        """
}

process prep_merge_samples_from_h5ad {
    input:
        tuple(
            val(experiment_id),
            path(file_h5ad),
        )
    label 'process_tiny'
    output:
        path("${experiment_id}---h5ad.h5ad", emit: h5ad)

    script:
        """
        ln --physical ${file_h5ad} ${experiment_id}---h5ad.h5ad
        """
}