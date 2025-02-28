
def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

process PCA {

    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
        // /software/hgi/containers/nf_qc_scrna_v1.img
    } else {
        container "mercury/nf_qc_scrna:v3"
    }

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("-", "")},
                mode: "${params.copy_mode}",
                overwrite: "true"
    
    publishDir  path: "${params.outdir}/handover/merged_h5ad/",
            saveAs: {filename ->
                    if (filename.contains("adata-normalized_pca-counts.h5ad")) {
                        filename = '5.adata-normalized_pca-counts.h5ad'
                    }else if (filename.contains("adata-normalized_pca-counts.h5ad")) {
                        filename = '5.adata-normalized_pca-counts.h5ad'
                    }  else {
                        null
                    }
                },
            mode: "${params.copy_mode}",
            overwrite: "true"  

    input:
        path(file__anndata)
        val(outdir)
        val(layer)

    output:
        val(outdir, emit: outdir)
        val("${outdir}", emit: outdir3)
        path("adata-normalized_pca.h5ad", emit: anndata)
        path("adata-metadata.tsv.gz", emit: metadata)
        path("adata-pcs.tsv.gz", emit: pcs)
        path(
            "adata-normalized_pca-counts.h5ad",
            emit: anndata_filtered_counts
        )
        val("${param_details}", emit: param_details)
        path("plots/*.pdf")
        path("plots/*.png") optional true

    script:

        """
        rm -fr plots
        pca_anndata.py \
            --h5_anndata ${file__anndata} \
            --overwrite_x_with_layer ${layer} \
            --output_file adata \
            --number_cpu ${task.cpus} \
            --drop_cell_passes_qc_from_clustering ${params.drop_cell_passes_qc_from_clustering}
        mkdir plots
        
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """

}


process NORMALISE_AND_PCA {
    // Takes annData object, nomalizes across samples, calculates PCs.
    // NOTE: Once normalization is set, it would be faster to normalize per
    //       sample and then merge.
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    // Couple of steps are happening here:
    // 1) Dropping any donors that do not pass the tresholds of number of cells.
    // 2) Regressing any variables that isers has specified.
    // 3) Filters out any unnecessary genes that users have specified.
    // 4) Scores the genes according to their gene score as per table provided.

    scratch false      // use tmp directory
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "mercury/nf_qc_scrna:v3"
    }


    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("-", "")},
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        
        path(file__anndata)
        val(analysis_mode)
        val(layer)
        val(minimum_number_of_cells_for_donor)
        val(file__genes_exclude_hvg)
        val(file__genes_score)
        val(genes_exclude)
        val(genes_at_least_in_nr_cells)
        each vars_to_regress

    output:
        
        val(outdir, emit: outdir)
        val("${outdir}", emit: outdir3)
        path("adata-normalized.h5ad", emit: anndata)
        val("${param_details}", emit: param_details)
        path("plots/*.pdf")
        path("plots/*.png") optional true
        path('donor_level_anndata_QCfiltered/*___sample_QCd_adata.h5ad',emit: sample_QCd_adata)

    script:
        
        analysis_mode = "${analysis_mode}"
        if (analysis_mode == "subclustering"){
            layer = "${layer}"
        }
        // Add any variables we are regressing to the output dir.
        param_details="vars_to_regress=none"
        if (vars_to_regress == "") {
            cmd__vars_to_regress = ""
        } else {
            param_details = "vars_to_regress=${vars_to_regress}"
            cmd__vars_to_regress = "--vars_to_regress ${vars_to_regress}"
        }

        outdir = "${params.outdir}/clustering_and_integration/normalize=total_count.${param_details}"
        // Add details on the genes we are exlcuding from hgv list.
        // file_vge = "${file__genes_exclude_hvg}"
        // outdir = "${outdir}.hvg_exclude=${file_vge}"
        // Add details on the scores we are using.
        // file_score = "${file__genes_score}"
        // outdir = "${outdir}.scores=${file_score}"
        // this is where the subfolder 1 is determined
        // Customize command for optional files.
        cmd__genes_exclude_hvg = ""
        if (file__genes_exclude_hvg != "") {
            cmd__genes_exclude_hvg = "--variable_genes_exclude '${file__genes_exclude_hvg}'"
        }

        cmd__genes_exclude = ""
        if (genes_exclude != "") {
            cmd__genes_exclude = "--exclude_gene_list '${genes_exclude}'"
        }
        
        cmd__genes_score = ""
        if (file__genes_score != "") {
            cmd__genes_score = "--score_genes ${file__genes_score}"
        }

        """
        rm -fr plots
        0035-scanpy_normalize_pca.py \
            --h5_anndata ${file__anndata} \
            --overwrite_x_with_layer ${layer} \
            --output_file adata \
            --number_cpu ${task.cpus} \
            --minimum_number_of_cells_for_donor ${minimum_number_of_cells_for_donor} \
            ${cmd__vars_to_regress} \
            ${cmd__genes_exclude_hvg} \
            ${cmd__genes_exclude} \
            ${cmd__genes_score} \
            --drop_cell_passes_qc_from_clustering ${params.normalise.drop_cell_passes_qc_from_clustering} \
            --genes_at_least_in_nr_cells ${genes_at_least_in_nr_cells}
        mkdir plots
        
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}
