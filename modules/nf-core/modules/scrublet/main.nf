

def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}

process SCRUBLET {
    // Runs scrublet for each sample.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run

    tag "${samplename}"
    scratch true        // use tmp directory
    echo echo_mode       // echo output from script

    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }
    
    publishDir  path: "${outdir}",
                saveAs: {filename ->
                    if (filename.endsWith("multiplet_calls_published.txt")) {
                        null
                    } else {
                        filename.replaceAll("${runid}-", "")
                    }
                },
                mode: "${params.copy_mode}",
                overwrite: "true"




    input:
        val(outdir_prev)
        tuple(
            val(experiment_id),
            path(file_10x_barcodes),
            path(file_10x_features),
            path(file_10x_matrix)
        )
        val(expected_multiplet_rate)
        val(n_simulated_multiplet)
        val(multiplet_threshold_method)
        val(scale_log10)

    output:
        val(outdir, emit: outdir)
        val(experiment_id, emit: experiment_id)
        path("${runid}-${outfile}-scrublet.tsv.gz", emit: multiplet_calls)
        path(
            "${runid}-${outfile}-multiplet_calls_published.txt",
            emit: multiplet_calls_published
        )
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}/multiplet"
        outdir = "${outdir}.method=scrublet"
        outfile = "${experiment_id}"
        // Check to see if we should use use log10 of the doublet simulations
        // to derive the threshold
        cmd__scale_log10 = ""
        if (scale_log10 == "True") {
            cmd__scale_log10 = "--scale_log10"
        }
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        echo "run_scrublet: ${process_info}"
        echo "publish_directory: ${outdir}"
        rm -fr plots
        TMP_DIR=\$(mktemp -d -p \$(pwd))
        ln --physical ${file_10x_barcodes} \$TMP_DIR
        ln --physical ${file_10x_features} \$TMP_DIR
        ln --physical ${file_10x_matrix} \$TMP_DIR
        run_scrublet.py \
            --tenxdata_dir \$TMP_DIR \
            --expected_multiplet_rate ${expected_multiplet_rate} \
            --n_simulated_multiplet ${n_simulated_multiplet} \
            --multiplet_threshold_method ${multiplet_threshold_method} \
            ${cmd__scale_log10} \
            --output_file ${runid}-${outfile}
        echo -e "${experiment_id}\t${outdir}/${outfile}-scrublet.tsv.gz" > \
            ${runid}-${outfile}-multiplet_calls_published.txt
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true
        """
}
