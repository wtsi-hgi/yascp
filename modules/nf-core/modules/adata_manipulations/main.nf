process ADD_EXTRA_METADATA_TO_H5AD{
    tag "${samplename}"    
    label 'process_low'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
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
