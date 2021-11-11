process SPLIT_BATCH_H5AD {
    tag "${samplename}"    
    label 'process_low'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }
    
    input:
        path(file__anndata) // anndata h5ad file seurat_azimuth_pbmc_1.0

    output:
        path("${outfil_prfx}_*.h5ad", emit:files_anndata_batch)
        path("Samples.tsv", emit:sample_file)
        path(outfile, emit: file_batch_list)

    script:
        process_info = "${task.cpus} (cpus), ${task.memory} (memory)"
        outfil_prfx = "${file__anndata}".minus(".h5ad")
        outfile = "${outfil_prfx}".plus("_files.txt")
        """
           scanpy_split_h5ad.py ${file__anndata} ${outfil_prfx}
        """

}