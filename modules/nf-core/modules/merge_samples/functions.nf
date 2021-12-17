def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

process dummy_filtered_channel{
    label 'process_low' 
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/wtsihgi_nf_scrna_qc_0417190-2021-12-16-133460e8fb0b.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
    } else {
        container "wtsihgi/nf_scrna_qc:0417190"
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
    //tag { output_dir }
    //cache true        // cache results from run
    tag "${samplename}"
    
    label 'process_low' 
    publishDir  path: "${outdir}/merged_h5ad",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${params.copy_mode}",
                overwrite: "true"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/wtsihgi_nf_scrna_qc_0417190-2021-12-16-133460e8fb0b.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
    } else {
        container "wtsihgi/nf_scrna_qc:0417190"
        //// container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }

    input:
        val(outdir_prev)
        path(file_paths_h5ad)
        path(file_metadata)
        val(file_params)
        val(file_cellmetadata)
        val(metadata_key)
        file(file_h5ad)
        val(anndata_compression_opts)

    // NOTE: use path here and not file see:
    //       https://github.com/nextflow-io/nextflow/issues/1414
    output:
        path("${runid}-adata.h5ad", emit: anndata)
        path(
            "${runid}-adata-cell_filtered_per_experiment.tsv.gz",
            emit: cells_filtered
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // String filename = './parameters.yml'
        // yaml.dump(file_params , new FileWriter(filename))
        // Customize command for optional files.
        cmd__params = ""
        if (file_params != "no_file__file_sample_qc") {
            cmd__params = "--params_yaml ${file_params}"
        }
        cmd__cellmetadata = ""
        if (file_cellmetadata != "no_file__file_cellmetadata") {
            cmd__cellmetadata = "--cell_metadata_file ${file_cellmetadata}"
        }
        files__h5ad = file_h5ad.join(',')
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "merge_samples: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        nf_helper__prep_h5addata_file.py \
            --h5ad_list ${files__h5ad} \
            --output_file nf_prepped__file_paths_h5ad.tsv
        scanpy_merge_from_h5ad.py \
            --h5addata_file nf_prepped__file_paths_h5ad.tsv \
            --sample_metadata_file ${file_metadata} \
            --sample_metadata_columns_delete "sample_status,study,study_id" \
            --metadata_key ${metadata_key} \
            --number_cpu ${task.cpus} \
            --output_file ${runid}-adata \
            --anndata_compression_opts ${anndata_compression_opts} \
            ${cmd__params} \
            ${cmd__cellmetadata}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}

process merge_samples {
    // Takes a list of raw 10x files and merges them into one anndata object.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache true        // cache results from run
    
    tag "${samplename}"
    
    label 'process_low'

    publishDir  path: "${outdir}/merged_h5ad",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${params.copy_mode}",
                overwrite: "true"
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/wtsihgi_nf_scrna_qc_0417190-2021-12-16-133460e8fb0b.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
    } else {
        container "wtsihgi/nf_scrna_qc:0417190"
        //// container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }

    input:
        val(outdir_prev)
        path(file_paths_10x)
        path(file_metadata)
        val(file_params)
        val(file_cellmetadata)
        val(metadata_key)
        file(file_10x_barcodes)
        file(file_10x_features)
        file(file_10x_matrix)
        val(anndata_compression_opts)

    // NOTE: use path here and not file see:
    //       https://github.com/nextflow-io/nextflow/issues/1414
    output:
        path("${runid}-adata.h5ad", emit: anndata)
        path(
            "${runid}-adata-cell_filtered_per_experiment.tsv.gz",
            emit: cells_filtered
        )
        path("plots/*.png") optional true
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // String filename = './parameters.yml'
        // yaml.dump(file_params , new FileWriter(filename))
        // Customize command for optional files.
        cmd__params = ""
        if (file_params != "no_file__file_sample_qc") {
            cmd__params = "--params_yaml ${file_params}"
        }
        cmd__cellmetadata = ""
        if (file_cellmetadata != "no_file__file_cellmetadata") {
            cmd__cellmetadata = "--cell_metadata_file ${file_cellmetadata}"
        }
        files__barcodes = file_10x_barcodes.join(',')
        files__features = file_10x_features.join(',')
        files__matrix = file_10x_matrix.join(',')
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "merge_samples: ${process_info}"
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
            --output_file ${runid}-adata \
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

    output:
        path("${experiment_id}---h5ad.h5ad", emit: h5ad)

    script:
        """
        ln --physical ${file_h5ad} ${experiment_id}---h5ad.h5ad
        """
}
