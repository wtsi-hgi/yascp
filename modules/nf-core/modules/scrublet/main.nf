

def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}

process SCRUBLET {
    // Runs scrublet for each sample.
    // ------------------------------------------------------------------------
    //cache false        // cache results from run

    tag "${experiment_id}"
    //scratch true        // use tmp directory


    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }
    
    publishDir  path: "${params.outdir}/doublets/multiplet.method=scrublet",
                saveAs: {filename ->
                    if (filename.endsWith("multiplet_calls_published.txt")) {
                        null
                    }
                    else if (filename.endsWith(".gz")){
                        null
                    }
                    else {
                        filename.replaceAll("-", "")
                    }
                },
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        // path(outdir_prev)
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
        // val(outdir, emit: outdir)
        tuple val(experiment_id), path("${outfile}-scrublet.tsv.gz"), emit: scrublet_paths
        val(experiment_id, emit: experiment_id)
        path("${outfile}-scrublet.tsv.gz", emit: multiplet_calls)
        path(
            "${outfile}-multiplet_calls_published.txt",
            emit: multiplet_calls_published
        )
        path("plots/*.pdf") optional true
        path("plots/*.png") optional true
        tuple val(experiment_id), path("${experiment_id}-scrublet.tsv"), emit: result

    script:
        
        outdir = "${params.outdir}/doublets/multiplet"
        outdir = "${outdir}.method=scrublet"
        outfile = "${experiment_id}"
        // Check to see if we should use use log10 of the doublet simulations
        // to derive the threshold
        cmd__scale_log10 = ""
        if (scale_log10 == "True") {
            cmd__scale_log10 = "--scale_log10"
        }

        """

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
            --output_file ${outfile}
        echo -e "${experiment_id}\t${outdir}/${outfile}-scrublet.tsv.gz" > \
            ${outfile}-multiplet_calls_published.txt
        mkdir plots
        mv *pdf plots/ 2>/dev/null || true
        mv *png plots/ 2>/dev/null || true

        zcat ${experiment_id}-scrublet.tsv.gz | awk -F' ' '{print \$1"\\t"\$3"\\t"\$2}' > ${experiment_id}-scrublet.tsv

        """
}
