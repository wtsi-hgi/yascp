process CELLTYPE_FILE_MERGE{
    tag "${samplename}"    
    label 'process_high'
    publishDir  path: "${params.outdir}/celltype_assignemt/",
            saveAs: {filename ->
                    if (filename.contains("adata.h5ad")) {
                        null
                    } else {
                        filename
                    }
                },
            mode: "${params.copy_mode}",
            overwrite: "true"  

    publishDir  path: "${params.outdir}/handover/merged_h5ad/",
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
        container "${params.yascp_container}"
    } else {
       container "${params.yascp_container_docker}"
    }
    output:
        // path('adata.h5ad', emit:file__anndata_merged2)
        path("All_Celltype_Assignments.tsv",emit:celltype_assignments)
        path "tranche_celltype_report.tsv"
        path "donor_celltype_report.tsv"

    input:
        tuple path(azimuth_files), path(celltypist_paths), path(all_other_paths)
        tuple val(expid), path(file__anndata_input)
    script:
        def merged_files_outpath = workflow.workDir.toString()
        file(merged_files_outpath).mkdirs()
        def azimuth_files_path = "${merged_files_outpath}/azimuth_files.tsv"
        def celltypist_files_path = "${merged_files_outpath}/celltypist_files.tsv"
        def all_other_files_path = "${merged_files_outpath}/other_files.tsv"
        def adatas_path = "${merged_files_outpath}/adatas.tsv"

        new File(azimuth_files_path).text = azimuth_files.join("\n")
        new File(celltypist_files_path).text = celltypist_paths.join("\n")

        if ("${all_other_paths}" != 'fake_file.fq') {
            new File(all_other_files_path).text = all_other_paths.join("\n")
            other_paths = "--all_other_paths ${all_other_files_path}"
        } else {
            other_paths = ""
        }

        new File(adatas_path).text = file__anndata_input.join("\n")

        """
        generate_combined_celltype_anotation_file.py --all_azimuth_files ${azimuth_files_path} --all_celltypist_files ${celltypist_files_path} ${other_paths} --adata '${adatas_path}'
        """

}
