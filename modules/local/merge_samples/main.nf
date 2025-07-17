include { merge_samples_from_h5ad; merge_samples;prep_merge_samples;prep_merge_samples_from_h5ad } from './functions'

params.file_cellmetadata    = "no_file__file_cellmetadata"
params.metadata_key_column = [
    value: "experiment_id"
]
params.anndata_compression_opts = 9

process HASTAG_FILE_MERGE{
    tag "${samplename}"    
    label 'process_low'

    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
       container "${params.yascp_container_docker}"
    }
    output:
        path("All_hastag_Assignments.tsv",emit:results)

    input:
        path(hastag_files)

    script:
        def merged_files_outpath = workflow.workDir.toString()
        file(merged_files_outpath).mkdirs()
        def hastag_files_path = "${merged_files_outpath}/hastag_files.tsv"
        new File(hastag_files_path).text = hastag_files.join("\n")

        """
            generate_combined_hastag_anotation_file.py --all_hastag_files ${hastag_files_path}
        """

}


workflow MERGE_SAMPLES{
    take:
        channel__file_paths_10x
        file_metadata
        celltype_file
        hastag_labels
        mode
    main:
        log.info """---Merging samples in a single h5ad file---"""

        if (mode == 'h5ad'){
            log.info """---Using h5ad input file outputs from deconvolution---"""
            channel__file_paths_10x.splitCsv(header: true, sep: "\t", by: 1)
                .map{row -> tuple(row.experiment_id, file(row.h5ad_filepath))}.set{channel__file_paths_10x_paths}
            prep_merge_samples_from_h5ad(channel__file_paths_10x_paths)
            input_merge =  prep_merge_samples_from_h5ad.out.h5ad.collect()

        }else if (mode == 'barcodes'){
            log.info """---Using barcode input---"""
            prep_merge_samples(channel__file_paths_10x)
            input_merge = prep_merge_samples.out.mtx.collect()

        }
        
        merge_samples_from_h5ad(
                file_metadata,
                params.file_cellmetadata,
                params.metadata_key_column.value,
                input_merge,
                celltype_file,
                hastag_labels
        )
        
        file__anndata_merged = merge_samples_from_h5ad.out.anndata

    emit:
        file__anndata_merged
}
