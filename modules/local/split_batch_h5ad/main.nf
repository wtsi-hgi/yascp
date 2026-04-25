process SPLIT_BATCH_H5AD {
    tag "${samplename}"    
    label 'process_low'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "${params.yascp_container_docker}"
    }
    
    input:
        path(file__anndata)
        val(doublet_celltype_split_column)

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
           scanpy_split_h5ad.py ${file__anndata} ${outfil_prfx} ${doublet_celltype_split_column}
        """

}


