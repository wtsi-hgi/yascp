
def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

process NORMALISE_AND_PCA {
    // Takes annData object, nomalizes across samples, calculates PCs.
    // NOTE: Once normalization is set, it would be faster to normalize per
    //       sample and then merge.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }
    

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path(file__anndata)
        val(analysis_mode)
        val(layer)
        path(file__genes_exclude_hvg)
        path(file__genes_score)
        each vars_to_regress

    output:
        
        val(outdir, emit: outdir)
        
        val("${outdir}", emit: outdir3)
        path("${runid}-adata-normalized_pca.h5ad", emit: anndata)
        path("${runid}-adata-metadata.tsv.gz", emit: metadata)
        path("${runid}-adata-pcs.tsv.gz", emit: pcs)
        
        path(
            "${runid}-adata-normalized_pca-counts.h5ad",
            emit: anndata_filtered_counts
        )
        val("${param_details}", emit: param_details)
        path("plots/*.pdf")
        path("plots/*.png") optional true
        // tuple(
        //     val(outdir),
        //     path("${runid}-adata-normalized_pca.h5ad"),
        //     path("${runid}-adata-metadata.tsv.gz"),
        //     path("${runid}-adata-pcs.tsv.gz"),
        //     emit: results
        // )

    script:
        runid = random_hex(16)
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

        // todo - mo11 - these paths are confusing

        outdir = "${outdir_prev}/normalize=total_count.${param_details}"
        // Add details on the genes we are exlcuding from hgv list.
        file_vge = "${file__genes_exclude_hvg.getSimpleName()}"
        outdir = "${outdir}.hvg_exclude=${file_vge}"
        // Add details on the scores we are using.
        file_score = "${file__genes_score.getSimpleName()}"
        outdir = "${outdir}.scores=${file_score}"


        // this is where the subfolder 1 is determined

        // Customize command for optional files.
        cmd__genes_exclude_hvg = ""
        if (file__genes_exclude_hvg.name != "no_file__genes_exclude_hvg") {
            cmd__genes_exclude_hvg = "--variable_genes_exclude ${file__genes_exclude_hvg}"
        }
        cmd__genes_score = ""
        if (file__genes_score.name != "no_file__genes_score") {
            cmd__genes_score = "--score_genes ${file__genes_score}"
        }
        // Basic details on the run.
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"

        """
        rm -fr plots
        0035-scanpy_normalize_pca.py \
            --h5_anndata ${file__anndata} \
            --overwrite_x_with_layer ${layer} \
            --output_file ${runid}-adata \
            --number_cpu ${task.cpus} \
            ${cmd__vars_to_regress} \
            ${cmd__genes_exclude_hvg} \
            ${cmd__genes_score}
        mkdir plots
        
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
        // Old version with bash evaluation of optional commands
        //
        // echo "normalize_pca: ${process_info}"
        // # If there are entries in the variable_genes_exclude file, add it to
        // # the call.
        // cmd__vg_exclude="--variable_genes_exclude ${file__genes_exclude_hvg}"
        // val=\$(cat ${file__genes_exclude_hvg} | wc -l)
        // if [ \$val -eq 0 ]; then cmd__vg_exclude=""; fi
        // # If there are entries in the score_genes file, add it to the call.
        // cmd__score_genes="--score_genes ${file__genes_score}"
        // val=\$(cat ${file__genes_score} | wc -l)
        // if [ \$val -eq 0 ]; then cmd__score_genes=""; fi
        // 0035-scanpy_normalize_pca.py \
        //     --h5_anndata ${file__anndata} \
        //     --output_file ${runid}-adata \
        //     --number_cpu ${task.cpus} \
        //     ${cmd__vars_to_regress} \
        //     \${cmd__vg_exclude} \
        //     \${cmd__score_genes}
        // mkdir plots
        // mv *pdf plots/ 2>/dev/null || true
        // mv *png plots/ 2>/dev/null || true
}