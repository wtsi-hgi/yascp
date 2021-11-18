process ADD_EXTRA_METADATA_TO_H5AD{
    tag "${samplename}"    
    label 'process_low'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }


    input:
        path(file__anndata) // anndata h5ad file seurat_azimuth_pbmc_1.0
        path(metadata)
    
    output:
        path("andata_with_metadata.h5ad"), emit: file__anndata

    script:
    """
        add_extra_metadata.py --andata ${file__anndata} --metadata ${metadata} --anndata_compression_level ${params.split_h5ad_per_donor.anndata_compression_level}
    """

}