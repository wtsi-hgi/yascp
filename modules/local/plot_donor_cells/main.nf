
process PLOT_DONOR_CELLS {

    tag "${sample_donor_summary_tsv}"
    
    label 'process_low'
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    publishDir "${params.outdir}/plots/", mode: "${params.plot_donor_ncells.copy_mode}", overwrite: true,
        saveAs: { filename -> 
            if (filename == 'versions.yaml') {
                null
            } else if (filename.indexOf(".pdf") > 0) {
                filename.replaceFirst("outputs/","")
            } else {
                filename
            }
        }
    
    when: 
    params.plot_donor_ncells.run

    input: 
    path(sample_donor_summary_tsv)

    output: 
    path("outputs/*.pdf"), emit: sample_pdf
    path "versions.yml", emit: versions

    script:
    """
        python plot_donor_ncells.py \\
        --output_dir \$PWD/outputs \\
        --sample_donor_summary_tsv ${sample_donor_summary_tsv} \\
        --plotnine_dpi ${params.plot_donor_ncells.plotnine_dpi}

        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            python library argparse: \$(python -c "import argparse; print(argparse.__version__)")
            python library csv: \$(python -c "import csv; print(csv.__version__)")
            python library click: \$(python -c "import click; print(click.__version__)")
            python library matplotlib: \$(python -c "import matplotlib; print(matplotlib.__version__)")
            python library logging: \$(python -c "import logging; print(logging.__version__)")
            python library numpy: \$(python -c "import numpy; print(numpy.__version__)")
            python library pandas: \$(python -c "import pandas; print(pandas.__version__)")
            python library plotnine: \$(python -c "import plotnine; print(plotnine.__version__)")
            python library scanpy: \$(python -c "import scanpy; print(scanpy.__version__)")
            python library seaborn: \$(python -c "import seaborn; print(seaborn.__version__)")
        END_VERSIONS
    """
}
