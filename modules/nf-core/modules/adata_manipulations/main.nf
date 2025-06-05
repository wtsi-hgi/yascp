process ADD_EXTRA_METADATA_TO_H5AD{
    tag "${samplename}"    
    label 'process_low'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }


    input:
        path(file_metadata) // anndata h5ad file seurat_azimuth_pbmc_1.0
        path(metadata)
        path(cellranger_metadata)
        path(sample_metadata)
    output:
        path("metadata_combined.csv"), emit: file__anndata

    script:
    """
        cp ${file_metadata} metadata_combined.csv
    """

}
