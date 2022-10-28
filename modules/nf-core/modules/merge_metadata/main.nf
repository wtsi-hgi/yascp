process prep_collectmetadata{
    label 'process_tiny'
    input:
        tuple val(experiment_id), path(metadata_path)
    output:
        path("${experiment_id}---metadata.csv", emit: metadata)
    script:
        """
            ln --physical ${metadata_path} ${experiment_id}---metadata.csv
        """
}

process merge_metadata{
    publishDir "${versionsDir}", pattern: "*.versions.yml", mode: "${params.versions.copy_mode}"

     input:
        file(file_metadata)
    output:
        path("full_metadata.tsv", emit: metadata)
        path ('*.versions.yml')         , emit: versions 

    script:
        files__metadata = file_metadata.join(',')
        """
            python ${projectDir}/bin/combine_metadata.py -d ${files__metadata}

            ####
    ## capture software version
    ####
    version=\$(python -V |sed "s/Python //g")
    echo "${task.process}:" > ${task.process}.versions.yml
    echo "    python: \$version" >> ${task.process}.versions.yml
        """
}
