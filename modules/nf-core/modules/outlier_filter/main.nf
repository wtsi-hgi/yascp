def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

process MERGE_OUTLIER_FILES{
    // Takes annData object, plots predicted outlier cells
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    tag "${samplename}"

    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"

    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    publishDir  path: "${params.outdir}/merged_h5ad/",
                saveAs: {filename ->
                    if (filename.contains("___sample_QCd_adata.h5ad")) {
                        null
                    } else if(filename.contains("outlier_filtered_adata.h5ad"))  {
                        filename = "4.outlier_filtered_adata.h5ad"
                    }else{
                        filename
                    }
                },
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        path(file__anndata)
        path(file__anndata_files_filtered)
        
    output:
        path('outlier_filtered_adata.h5ad', emit: anndata)
    script:
        """
            merge_outliers.py -h5 ${file__anndata}
        """    
    
}

process OUTLIER_FILTER {
    // Takes annData object, plots predicted outlier cells
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    tag "${samplename}"

    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"

    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    publishDir  path: "${outdir}/merged_h5ad/",
                saveAs: {filename ->
                    if (filename.contains("___sample_QCd_adata.h5ad")) {
                        null
                    } else if(filename.contains("outlier_filtered_adata.h5ad"))  {
                        filename = "4.outlier_filtered_adata.h5ad"
                    }else{
                        null
                    }
                },
                mode: "${params.copy_mode}",
                overwrite: "true"



    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__cells_filtered)
        val(outliers_fraction)
        val(max_samples)
        val(anndata_compression_opts)
        path(gt_outlier_input)
        val(gt_match_based_adaptive_qc_exclusion_pattern)
        each refset

    output:
        // path("merged_h5ad/outlier_filtered_adata.h5ad", emit: anndata)
        // path('donor_level_anndata_QCfiltered/*___sample_QCd_adata.h5ad',emit: sample_QCd_adata)
        path(
            "outlier_filtered_adata-outliers_filtered*.tsv",
            emit: cells_filtered
        )
        path("plots/*")
        // path("merged_h5ad/*")


    script:
        if(gt_match_based_adaptive_qc_exclusion_pattern!=''){
            filter_strategy_exclusion = "--patterns_exclude='${gt_match_based_adaptive_qc_exclusion_pattern}' --gt_match_file ${gt_outlier_input}"
        }else{
            filter_strategy_exclusion = ""
        }
        outdir = "${outdir_prev}"
        // Append run_id to output file.
        outfile = "outlier_filtered_adata"
        """
            echo "publish_directory: ${outdir}"
            rm -fr plots
            0026-filter_outlier_cells.py \
                --h5_anndata ${file__anndata} \
                --cell_filtered_per_experiment_file ${file__cells_filtered} \
                --outliers_fraction 0 \
                --metadata_columns ${refset.metadata_columns} \
                --cell_qc_column cell_passes_qc \
                --method ${refset.method} \
                --outliers_fraction ${outliers_fraction} \
                --max_samples ${max_samples} \
                --output_file ${outfile} \
                --anndata_compression_opts ${anndata_compression_opts} \
                --filter_strategy='${refset.outlier_filtering_strategy}' \
                --MAD_thresholds='${refset.mad_tresholds}' \
                ${filter_strategy_exclusion}
                
            mkdir plots
            mv *pdf plots/ 2>/dev/null || true
            mv *png plots/ 2>/dev/null || true
            mv per_celltype_outliers* plots/ 2>/dev/null || true
        """
}
