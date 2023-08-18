process prep_collectmetadata{
    label 'process_tiny'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/nf_scrna_deconv_v2.img"
        //// container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_latest.img"
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
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/nf_scrna_deconv_v2.img"
        //// container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_latest.img"
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
