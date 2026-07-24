process YASCP_INPUTS {
    // Takes annData object with PCs and returns plots
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    scratch false      // use tmp directory

    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }

    input:
        path(input_file)


    // the output for this is a correct format input files as per cb 6.1
    output:
        path "input_file_corectly_formatted.tsv", emit: input_file_corectly_formatted
        path "versions.yml", emit: versions


    script:
        
        """
        link_files.py --input_table ${input_file} \
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            python: \$(python --version | sed 's/Python //g')
            python library argparse: \$(python -c "import argparse; print(argparse.__version__)")
            python library pandas: \$(python -c "import pandas; print(pandas.__version__)")
        END_VERSIONS
        """
}