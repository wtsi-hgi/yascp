#!/usr/bin/env nextflow

process sccaf_assess_clustering {
  
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
        
    } else {
        container "wtsihgi/nf_scrna_qc_scaaf:1c77f49"
    }


  publishDir path: "${outdir}",
             mode: "${params.copy_mode}",
             overwrite: "true"

  when:
    params.sccaf.run_assessment

  input:
    val(outdir_prev)
    path(file__anndata) // anndata h5ad file
    path(file__external_clustering_tsv)

  output:
    path("${outfile_roc_pdf}")
    path("${outfile_prc_pdf}")
    path("${outfile_acc_txt}")

  script:
    outdir = "${outdir_prev}/sccaf/clustering_assessment"
    outfil_prfx = "${file__external_clustering_tsv}".minus(".tsv.gz").plus("_assess")
    outfile_roc_pdf = outfil_prfx.plus("_clust_roc.pdf")
    outfile_prc_pdf = outfil_prfx.plus("_clust_prc.pdf")
    outfile_acc_txt = outfil_prfx.plus("_acc.txt")
    """

    sccaf_assess_clustering.py \\
       --nthreads ${task.cpus} \\
       --use-pca \\
       --output-prefix ${outfil_prfx} \\
       ${file__anndata} \\
       ${file__external_clustering_tsv}
    """
}

process sccaf_optimize_clustering {


    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "wtsihgi/nf_scrna_qc_scaaf:1c77f49"
    }

  publishDir path: "${outdir}",
             saveAs: {filename -> outfil_prfx.plus(filename)},
             mode: "${params.copy_mode}",
             overwrite: "true"
             //saveAs: {filename -> outfil_prfx.plus(filename)},

  when:
    params.sccaf.run_optimization

  input:
    val(outdir_prev)
    path(file__anndata)
    path(file__external_clustering_tsv)
    val(min_accuracy)

    output:
      path("${outfile}")
      path("roc-curve.png")
      path("optim.png")
      path("rounds.txt")

  script:
    outdir = "${outdir_prev}/sccaf/clustering_optimization"
    outdir_rel = "optim_out"
    outfil_prfx = "${file__external_clustering_tsv}".minus(".tsv.gz").plus("_optim_")
    outfile = outfil_prfx.plus("clustopt.h5ad")
    """
        sccaf \\
            --output-file ${outfile} \\
            --cores ${task.cpus} \\
            --input-file ${file__anndata} \\
            --external-clustering-tsv ${file__external_clustering_tsv} \\
            --use-pca \\
            --optimise \\
            --min-accuracy ${min_accuracy} \\
            --produce-rounds-summary
        #      --optimisation-plots-output ${outdir_rel}
    """
}
