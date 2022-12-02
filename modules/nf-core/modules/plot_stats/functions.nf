def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

process plot_filtered_cells {
    // Takes annData object, plots filtered cells
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    publishDir  path: "${outdir}/plots",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir)
        path(file__filtered_cells)

    output:
        path("*.png") optional true
        path("*.pdf") optional true

    script:
        
        process_info = ""
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_filtered_cells: ${process_info}"
        echo "publish_directory: ${outdir}"
        0026-plot_filtered_cells.py \
            --tsv_file ${file__filtered_cells} \
            --output_file adata-cell_filtered_per_experiment

        """
}



process plot_pcs {
    // Takes annData object with PCs and returns plots
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory

    label 'process_medium_memory'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    publishDir  path: "${outdir}/plots",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir)
        path(file__anndata)
        each n_pcs
        val(colors_quantitative)
        val(colors_categorical)

    output:
        val(outdir, emit: outdir)
        path("*.png"), emit: out_png
        path("*.pdf") optional true

    script:
        
        outfile = "pca"
        cmd__colors_quant = ""
        if (colors_quantitative != "") {
            cmd__colors_quant = "--colors_quantitative ${colors_quantitative}"
        }
        cmd__colors_cat = ""
        if (colors_categorical != "") {
            cmd__colors_cat = "--colors_categorical ${colors_categorical}"
        }
        // drop_legend_n = "-1"
        // if (cmd__colors_cat.contains("experiment_id")) {
        //     drop_legend_n = "8"
        // }
        process_info = ""
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "pca_plot: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        pca_plot.py \
            --h5_anndata ${file__anndata} \
            --num_pcs ${n_pcs} \
            ${cmd__colors_quant} \
            ${cmd__colors_cat} \
            --output_file ${outfile}
        """
}


process plot_predicted_sex {
    // Takes annData object, plots the predicted sex fron gene expression
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    publishDir  path: "${outdir}/plots",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir)
        path(file__anndata)

    output:
        path("*.png") optional true
        path("*.pdf") optional true

    script:

        outfile = "scatterplot-sex_sample_swap_check"
        process_info = ""
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_predicted_sex: ${process_info}"
        echo "publish_directory: ${outdir}"
        0028-plot_predicted_sex.py \
            --h5_anndata ${file__anndata} \
            --output_file ${outfile}
        """
}

process plot_qc {
    // Takes annData object, generates basic qc plots
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
	
        container "/software/hgi/containers/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// working but deprecated: container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
	
	//// container "wtsihgi/nf_scrna_qc:6bb6af5"
        //// 0417190 fails:
	////   with error:     File "/home/container_user/conda/envs/cenv/lib/python3.9/site-packages/numpy/core/shape_base.py", line 283, in vstack               
            // return _nx.concatenate(arrs, 0)                                    
            // File "<__array_function__ internals>", line 5, in concatenate     
            // ValueError: all the input arrays must have same number of dimensions, but the array at index 0 has 2 dimension(s) and the array at index 2 has 3 dimension(s)
	
    }

    publishDir  path: "${outdir}/plots",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir)
        path(file__anndata)
        each facet_columns

    output:
        path("*.png")
        path("*.pdf") optional true
        path("*.tsv") optional true

    script:

        outfile = "outfile"
        // Figure out if we are facetting the plot and update accordingly.
        cmd__facet_columns = ""
        if (facet_columns != "") {
            cmd__facet_columns = "--facet_columns ${facet_columns}"
        }
        process_info = ""
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
            echo "plot_qc: ${process_info}"
            echo "publish_directory: ${outdir}"
            plot_qc_umi_nfeature_mt.py \
                --h5_anndata ${file__anndata} \
                --output_file ${outfile} \
                ${cmd__facet_columns}
            plot_qc_umi_mt_density.py \
                --h5_anndata ${file__anndata} \
                --output_file ${outfile} \
                ${cmd__facet_columns}
            plot_qc_nfeature_mt_density.py \
                --h5_anndata ${file__anndata} \
                --output_file ${outfile} \
                ${cmd__facet_columns}
            0027-calculate_mads.py \
                --h5_anndata ${file__anndata} \
                --qc_key 'pct_counts_gene_group__mito_transcript,pct_counts_gene_group__mito_protein,pct_counts_gene_group__ribo_protein,pct_counts_gene_group__ribo_rna,total_counts,n_genes_by_counts' \
                --output_file mads
        """
}

process plot_distributions {
    // Takes annData object, generates basic qc plots
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    tag "${samplename}"
    scratch true        // use tmp directory

    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    publishDir  path: "${outdir}/plots",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir)
        path(file__anndata)
        each facet_columns
        each variable_columns_distribution_plots

    output:

        path("*.png")
        path("*.pdf") optional true
        path("*.tsv") optional true

    script:
        outfile = "outfile"
        // Figure out if we are facetting the plot and update accordingly.
        cmd__facet_columns = ""
        if (facet_columns != "") {
            cmd__facet_columns = "--facet_columns ${facet_columns}"
        }
        // Run distribution across cells if a value is specified
        cmd__anndataobs = ""
        cmd__anndataobs_ecdf = ""
        if (variable_columns_distribution_plots != "") {
            cmd__anndataobs = "plot_anndataobs_distribution_across_cells.py"
            cmd__anndataobs = "${cmd__anndataobs} --h5_anndata ${file__anndata}"
            cmd__anndataobs = "${cmd__anndataobs} --output_file ${outfile}"
            cmd__anndataobs = "${cmd__anndataobs} --variable_columns ${variable_columns_distribution_plots}"
            cmd__anndataobs = "${cmd__anndataobs} ${cmd__facet_columns}"
            cmd__anndataobs_ecdf = "${cmd__anndataobs} --ecdf"
        }
        process_info = ""
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "plot_distributions: ${process_info}"
        echo "publish_directory: ${outdir}"
        ${cmd__anndataobs}
        ${cmd__anndataobs_ecdf}
        """
}
