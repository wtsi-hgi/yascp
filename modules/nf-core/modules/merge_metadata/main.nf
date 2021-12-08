process prep_collectmetadata{
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
     input:
        file(file_metadata)
    output:
        path("full_metadata.tsv", emit: metadata)
    script:
        files__metadata = file_metadata.join(',')
        """
            python ${projectDir}/bin/combine_metadata.py -d ${files__metadata}
        """
}