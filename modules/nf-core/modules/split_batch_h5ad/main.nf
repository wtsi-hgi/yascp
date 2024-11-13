process SPLIT_BATCH_H5AD {
    tag "${samplename}"    
    label 'process_low'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.nf_scrna_qc_sif_container}"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }
    
    input:
        tuple val(name), path(file__anndata) // anndata h5ad file seurat_azimuth_pbmc_1.0
        val(mode)

    output:
        path("AZ_${outfil_prfx}_*.h5ad", emit:files_anndata_batch)
        path("Samples.tsv", emit:sample_file)
        path("AZ_Samples.tsv", emit:az_sample_file)
        path(outfile, emit: file_batch_list)
        path(file__anndata,emit:adata)
        path("${outfil_prfx}_*.h5ad", emit:keras_outfile)
        

    script:
        outfil_prfx = "${file__anndata}".minus(".h5ad")
        outfile = "${outfil_prfx}".plus("_files.txt")
        """

           scanpy_split_h5ad.py ${file__anndata} ${outfil_prfx} ${mode}
        """

}


