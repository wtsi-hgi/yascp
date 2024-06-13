process CELLTYPE_FILE_MERGE{
    tag "${samplename}"    
    label 'process_medium'
    publishDir  path: "${params.outdir}/celltype/",
            saveAs: {filename ->
                    if (filename.contains("adata.h5ad")) {
                        null
                    } else {
                        filename
                    }
                },
            mode: "${params.copy_mode}",
            overwrite: "true"  

    publishDir  path: "${params.outdir}/merged_h5ad/",
            saveAs: {filename ->
                    if (filename.contains("adata.h5ad")) {
                        filename = "2.celltype_anotated_merged.h5ad"
                    } else {
                        null
                    }
                },
            mode: "${params.copy_mode}",
            overwrite: "true"  


    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
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
        path(azimuth_files)
        path(celltypist_paths)
        path(all_other_paths)
        path(file__anndata_input)
    script:
        all_azimuth_files = azimuth_files.join("::")
        all_celltypist_files = celltypist_paths.join("::")
        if ("${all_other_paths}"!='fake_file.fq'){
            all_other_paths_comb = all_other_paths.join("::")
            other_paths ="--all_other_paths ${all_other_paths_comb}"
        }else{
            other_paths = ""
        }
        
        all_adatas = file__anndata_input.join("::")
        """
            generate_combined_celltype_anotation_file.py --all_azimuth_files ${all_azimuth_files} --all_celltypist_files ${all_celltypist_files} ${other_paths} --adata '${all_adatas}'
        """

}
