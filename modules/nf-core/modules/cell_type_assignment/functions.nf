process CELLTYPE_FILE_MERGE{
    tag "${samplename}"    
    label 'process_low'
    publishDir  path: "${params.outdir}/celltype/",
            saveAs: {filename -> filename},
            mode: "${params.copy_mode}",
            overwrite: "true"  
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        // container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }
    output:
        path('adata.h5ad', emit:file__anndata_merged2)
        path("All_Celltype_Assignments.csv",emit:celltype_assignments)
        path "tranche_celltype_report.tsv"
        path "donor_celltype_report.tsv"

    input:
        file(azimuth_files)
        file(celltypist_paths)
        file(all_other_paths)
        path(file__anndata_input)
    script:
        all_azimuth_files = azimuth_files.join("::")
        all_celltypist_files = celltypist_paths.join("::")
        all_other_paths_comb = all_other_paths.join("::")
        """
            generate_combined_celltype_anotation_file.py --all_azimuth_files ${all_azimuth_files} --all_celltypist_files ${all_celltypist_files} --all_other_paths ${all_other_paths_comb} --adata ${file__anndata_input}
        """

}
