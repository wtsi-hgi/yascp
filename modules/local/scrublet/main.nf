

def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

if (binding.hasVariable("echo_mode") == false) {
    echo_mode = true
}

process SCRUBLET {
    // Runs scrublet for each sample.

    tag "${experiment_id}"

    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }
    
    publishDir  path: "${params.outdir}/doublet_detection/scrublet",
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
        path "versions.yml", emit: versions

    script:
        
        outdir = "${params.outdir}/doublet_detection/"
        outdir = "${outdir}scrublet"
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

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            scrublet: \$(python -c "import scrublet; print(scrublet.__version__)")
            scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
            argparse: \$(python -c "import argparse; print(argparse.__version__)")
            os: \$(python -c "import os; print(os.__version__)")
            csv: \$(python -c "import csv; print(csv.__version__)")
            random: \$(python -c "import random; print(random.__version__)")
            numpy: \$(python -c "import numpy; print(numpy.__version__)")
            pandas: \$(python -c "import pandas; print(pandas.__version__)")
            matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
            seaborn: \$(python -c "import seaborn; print(seaborn.__version__)")
            plotnine: \$(python -c "import plotnine; print(plotnine.__version__)")
            skimage: \$(python -c "import skimage; print(skimage.__version__)")
            distutils: \$(python -c "import distutils; print(distutils.__version__)")
        END_VERSIONS
        """
}
