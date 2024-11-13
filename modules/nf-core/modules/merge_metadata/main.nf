process prep_collectmetadata{
    label 'process_tiny'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.nf_scrna_qc_v3_container}"
    } else {
        container "mercury/scrna_deconvolution:v2"
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

process merge_metadata{
    label 'process_tiny'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.nf_scrna_qc_v3_container}"
    } else {
        container "mercury/scrna_deconvolution:v2"
    }

    input:
        file(file_metadata)
    output:
        path("full_metadata.tsv", emit: metadata)
    script:
        files__metadata = file_metadata.join(',')
        """
            combine_metadata.py -d ${files__metadata}
        """
}
