include { merge_samples_from_h5ad; merge_samples;prep_merge_samples;prep_merge_samples_from_h5ad } from './functions'

params.file_cellmetadata    = "no_file__file_cellmetadata"
params.metadata_key_column = [
    value: "experiment_id"
]
params.anndata_compression_opts = 9

workflow MERGE_SAMPLES{
    take:
        channel__file_paths_10x
        file_metadata
        celltype_file
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

        input_merge.subscribe { println "input_merge: $it" }
        channel__file_paths_10x.subscribe { println "channel__file_paths_10x: $it" }
        file_metadata.subscribe { println "file_metadata: $it" }
        celltype_file.subscribe { println "celltype_file: $it" }
        
        merge_samples_from_h5ad(
                file_metadata,
                params.file_cellmetadata,
                params.metadata_key_column.value,
                input_merge,
                celltype_file
        )
        
        file__anndata_merged = merge_samples_from_h5ad.out.anndata

    emit:
        file__anndata_merged
}
