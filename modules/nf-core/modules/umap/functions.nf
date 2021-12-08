#!/usr/bin/env nextflow

def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}

process umap_calculate {
    // UMAP from reduced_dims.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory

    // Don't publish these results as they are just temporary
    // publishDir  path: "${outdir}",
    //             saveAs: {filename -> filename.replaceAll("${runid}-", "")},
    //             mode: "${task.publish_mode}",
    //             overwrite: "true"
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }

    input:
        val(outdir_prev)
        path(file__anndata)
        path(file__metadata) // not really needed
        path(file__pcs) // not really needed
        path(file__reduced_dims)
        val(use_pcs_as_reduced_dims) // For BBKNN
        each n_neighbors
        each umap_init
        each umap_min_dist
        each umap_spread

    output:
        val(outdir, emit: outdir)
        // tuple(
        //     val("${in_file_id}"),
        //     file("${runid}-${outfile}.h5ad"),
        //     emit: anndata
        // )
        tuple(
            val("${in_file_id}"),
            val("${outdir}"),
            file(file__anndata),
            file(file__metadata),
            file(file__pcs),
            file(file__reduced_dims),
            file("${runid}-${outfile}.h5ad"),
            emit: outdir_anndata
        )

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // Get a unique identifier for the input file that will be later
        // used as a key to merge umap jobs on the same input.
        //
        // Here the unique identifier is file__reduced_dims, since the reduced
        // dims function does not change anndata but just writes a tsv file.
        //
        // If file__reduced_dims == file object, then use
        // file__reduced_dims.name...
        in_file_id = file__reduced_dims.toString().tokenize('-').get(0)
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad").split("-").drop(1).join("-")
        outfile = "${outfile}-umap"
        // Check to see if we should use PCs in the reduced dims slot.
        // Important for BBKNN where reduced_dims == UMAPs not adjusted PCs.
        cmd__tsv_pcs = "--tsv_pcs ${file__reduced_dims}"
        if (use_pcs_as_reduced_dims == "True") {
            cmd__tsv_pcs = "--tsv_pcs ${file__pcs}"
        }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
            umap_calculate.py \
                --h5_anndata ${file__anndata} \
                ${cmd__tsv_pcs} \
                --n_neighbors ${n_neighbors} \
                --umap_init ${umap_init} \
                --umap_min_dist ${umap_min_dist} \
                --umap_spread ${umap_spread} \
                --number_cpu ${task.cpus} \
                --output_file ${runid}-${outfile}
        """
        //--calculate_densities \
}

process generate_final_UMAPS{

    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }

  publishDir  path: "${outdir}",
              mode: "${params.copy_mode}",
              overwrite: "true"
  input:
    path(file__anndata)
    val(outdir_prev)
  output:
    path("umap-*")
  
  script:
    outdir = "${outdir_prev}/UMAPs"
    """
        umap_plot_final.py \
            --h5_anndata ${file__anndata} \
            --number_cpu 1 \
            --colors_quantitative n_cells,total_counts,pct_counts_gene_group__mito_transcript,prob_doublet,pct_counts_gene_group__ribo_rna,Azimuth:predicted.celltype.l2.score,Azimuth:mapping.score \
            --colors_categorical experiment_id,Azimuth:predicted.celltype.l2,Celltypist:Immune_All_Low,Celltypist:Immune_All_High,Celltypist:Immune_Blood_High,Celltypist:Immune_Blood_Low,donor_id  \
            --drop_legend_n 40 \
            --output_file UMAP
    """

}

process umap_gather {
    // Merge UMAP from reduced_dims (reduce or gather).
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
                saveAs: {filename ->
                    if(filename.endsWith("metadata.tsv.gz")) {
                        null
                    } else if(filename.endsWith("pcs.tsv.gz")) {
                        null
                    } else if(filename.endsWith("reduced_dims.tsv.gz")) {
                        null
                    } else {
                        filename.replaceAll("${runid}-", "")
                    }
                },
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        // val(outdir_prev)
        // path(original__file__anndata)
        // path(original__file__metadata)
        // path(original__file__pcs)
        // path(original__file__reduced_dims)
        // //tuple(val(key), path(files__anndata))
        // tuple(val(key), val(outdir_prev_tuple), path(files__anndata))
        tuple(
            val(key),
            val(outdir_prev),
            path(original__file__anndata),
            path(original__file__metadata),
            path(original__file__pcs),
            path(original__file__reduced_dims),
            path(files__anndata)
        )
    output:
        val(outdir_prev, emit: outdir)
        path("${runid}-${outfile}.h5ad", emit: anndata)
        path(original__file__metadata, emit: metadata)
        path(original__file__pcs, emit: pcs)
        path(original__file__reduced_dims, emit: reduced_dims)

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"  // For some reason dir here messed up?
        //outdir_prev_tuple = outdir_prev_tuple.unique().join("")
        //outdir = "${outdir_prev_tuple}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${original__file__anndata}".minus(".h5ad").split("-")
            .drop(1).join("-")
        outfile = "${outfile}-umap"
        // outfile = "adata-umap"
        // Get all of the adata files that we want to gather
        files__anndata = files__anndata.join(',')
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        // Two ways to call --h5_anndata_list
        // --h5_anndata_list ${files__anndata}
        // --h5_anndata_list \$(ls -1 *h5ad | paste -sd "," -)
        //sleep 5m
        """
       
        umap_gather.py \
            --h5_root ${original__file__anndata} \
            --output_file ${runid}-${outfile} \
            --h5_anndata_list ${files__anndata}
        """
}


process umap_plot_swarm {
    // Plot UMAPs.
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
        val(colors_quantitative)
        val(colors_categorical)
        val(drop_legend_n)

    output:
        path("plots/*.png")
        path("plots/*.pdf") optional true
        // path("plots/*.svg") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        // outfile = "${file__anndata}".minus(".h5ad").split("-").drop(1).join("-")
        outfile = "umap"
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
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        rm -fr plots
        umap_plot.py \
            --h5_anndata ${file__anndata} \
            --number_cpu ${task.cpus} \
            ${cmd__colors_quant} \
            ${cmd__colors_cat} \
            --drop_legend_n ${drop_legend_n} \
            --output_file ${runid}-${outfile}

        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}


process umap_calculate_and_plot {
    // UMAP from reduced_dims.
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
        path(file__pcs)
        path(file__reduced_dims)
        val(use_pcs_as_reduced_dims) // For BBKNN
        val(colors_quantitative)
        val(colors_categorical)
        each n_neighbors
        each umap_init
        each umap_min_dist
        each umap_spread

    output:
        path("plots/*.png")
        path("plots/*.pdf") optional true
        // path("plots/*.svg") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        // outfile = "${file__anndata}".minus(".h5ad").split("-").drop(1).join("-")
        outfile = "umap"
        cmd__colors_quant = ""
        if (colors_quantitative != "") {
            cmd__colors_quant = "--colors_quantitative ${colors_quantitative}"
        }
        cmd__colors_cat = ""
        if (colors_categorical != "") {
            cmd__colors_cat = "--colors_categorical ${colors_categorical}"
        }
        drop_legend_n = "-1"
        if (cmd__colors_cat.contains("experiment_id")) {
            drop_legend_n = "8"
        }
        // Check to see if we should use PCs in the reduced dims slot.
        // Important for BBKNN where reduced_dims == UMAPs not adjusted PCs.
        cmd__tsv_pcs = "--tsv_pcs ${file__reduced_dims}"
        if (use_pcs_as_reduced_dims == "True") {
            cmd__tsv_pcs = "--tsv_pcs ${file__pcs}"
        }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        rm -fr plots
        umap_calculate_and_plot.py \
            --h5_anndata ${file__anndata} \
            ${cmd__tsv_pcs} \
            --n_neighbors ${n_neighbors} \
            --umap_init ${umap_init} \
            --umap_min_dist ${umap_min_dist} \
            --umap_spread ${umap_spread} \
            --number_cpu ${task.cpus} \
            ${cmd__colors_quant} \
            ${cmd__colors_cat} \
            --drop_legend_n ${drop_legend_n} \
            --output_file ${runid}-${outfile}
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
        // --calculate_densities
}
