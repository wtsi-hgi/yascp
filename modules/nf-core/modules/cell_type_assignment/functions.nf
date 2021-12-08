process CELLTYPE_FILE_MERGE{
    tag "${samplename}"    
    label 'process_low'
    publishDir  path: "${params.outdir}/celltype/",
            saveAs: {filename -> filename},
            mode: "${params.copy_mode}",
            overwrite: "true"  
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "quay.io/biocontainers/multiqc:1.10.1--py_0"
    }
    output:
        path('adata.h5ad', emit:file__anndata_merged2)
        path "All_Celltype_Assignments.csv"

    input:
        file(azimuth_files)
        file(celltypist_paths)
        file(file__anndata_merged)
    script:
        all_azimuth_files = azimuth_files.join("::")
        all_celltypist_files = celltypist_paths.join("::")
        """
            generate_combined_celltype_anotation_file.py --all_azimuth_files ${all_azimuth_files} --all_celltypist_files ${all_celltypist_files} --adata ${file__anndata_merged}
        """

}