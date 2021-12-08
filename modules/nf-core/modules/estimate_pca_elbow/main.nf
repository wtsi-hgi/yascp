
def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}

process ESTIMATE_PCA_ELBOW {
    // Takes annData object, estiamtes the elbow in PC var explained.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run
    scratch false      // use tmp directory
    echo echo_mode          // echo output from script

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
        val(add_n_to_estimate)

    output:
        val(outdir, emit: outdir)
        path("${runid}-${outfile}.tsv", emit: pca_elbow_estimate)
        env(AUTO_ELBOW, emit: auto_elbow)
        path("plots/*.png")
        path("plots/*.pdf") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        log.info("""outdir = ${outdir}""")
        // For output file, use anndata name. First need to drop the runid
        // from the file__anndata job.
        outfile = "${file__anndata}".minus(".h5ad")
            .split("-").drop(1).join("-")
        outfile = "${outfile}-knee"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
            rm -fr plots
            0030-estimate_pca_elbow.py \
                --h5_anndata ${file__anndata} \
                --add_n_pcs_to_elbow ${add_n_to_estimate} \
                --output_file ${runid}-${outfile}
            mkdir plots
            mv *pdf plots/ 2>/dev/null || true
            mv *png plots/ 2>/dev/null || true
            AUTO_ELBOW=\$(cat ${runid}-${outfile}-auto_elbow_estimate.tsv)
        """
}