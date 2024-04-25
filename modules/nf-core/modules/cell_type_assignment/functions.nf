process CELLTYPE_FILE_MERGE{
    tag "${samplename}"    
    label 'process_medium'
    publishDir  path: "${params.outdir}/celltype/",
            saveAs: {filename -> filename},
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
        def azimuth_files_path = "${workDir}/azimuth_files.tsv"
        azimuth_files.each { file -> file.withWriterAppend { w -> w.println("${file}") } }

        def celltypist_files_path = "${workDir}/celltypist_files.tsv"
        celltypist_paths.each { file -> file.withWriterAppend { w -> w.println("${file}") } }

        def all_other_files_path = "${workDir}/other_files.tsv"
        if (all_other_paths != "fake_file.fq") {
            all_other_paths.each { file -> file.withWriterAppend { w -> w.println("${file}") } }
            other_paths_option = "--all_other_paths ${all_other_files_path}"
        } else {
            other_paths_option = ""
        }

        def adatas_path = "${workDir}/adatas.tsv"
        file__anndata_input.each { file -> file.withWriterAppend { w -> w.println("${file}") } }

        """
        generate_combined_celltype_annotation_file.py --all_azimuth_files ${azimuth_files_path} --all_celltypist_files ${celltypist_files_path} ${other_paths_option} --adata '${adatas_path}'
        """

}
