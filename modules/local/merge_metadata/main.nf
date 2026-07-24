process PREP_COLLECTMETADATA{
    label 'process_tiny'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }
    
    input:
        tuple val(experiment_id), path(metadata_path)
    output:
        path("${experiment_id}---metadata.csv", emit: metadata)
    script:
        """
            ln -s ${metadata_path} ${experiment_id}---metadata.csv
        """
}

process MERGE_METADATA{
    label 'process_tiny'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
       container "${params.yascp_container_docker}"
    }

    input:
        file(file_metadata)
    output:
        path "full_metadata.tsv", emit: metadata
        path "versions.yml", emit: versions
    script:
        files__metadata = file_metadata.join(',')
        """
            combine_metadata.py -d ${files__metadata}

            cat <<-END_VERSIONS > versions.yml
            "${task.process}":
                python: \$(python --version | sed 's/Python //g')
                python library argparse: \$(python -c "import argparse; print(argparse.__version__)")
                python library pandas: \$(python -c "import pandas; print(pandas.__version__)")
            END_VERSIONS
        """
}
