include { merge_samples_from_h5ad; merge_samples;prep_merge_samples;prep_merge_samples_from_h5ad } from './functions'

params.file_sample_qc       = "no_file__file_sample_qc"
params.file_cellmetadata    = "no_file__file_cellmetadata"
params.metadata_key_column = [
    value: "experiment_id"
]
params.anndata_compression_opts = 9

workflow MERGE_SAMPLES{
    take:
        channel__file_paths_10x
        file_metadata
        mode
    main:
        log.info """---Merging samples in a single h5ad file---"""

        if (mode == 'h5ad'){
            log.info """---Using h5ad input file outputs from deconvolution---"""
            
            channel__file_paths_10x.splitCsv(header: true, sep: "\t", by: 1)
                .map{row -> tuple(row.experiment_id, file(row.h5ad_filepath))}.set{channel__file_paths_10x_paths}
            
            file_metadata.view()
            prep_merge_samples_from_h5ad(channel__file_paths_10x_paths)
            merge_samples_from_h5ad(
                    params.output_dir,
                    channel__file_paths_10x,
                    file_metadata,
                    params.file_sample_qc,
                    params.file_cellmetadata,
                    params.metadata_key_column.value,
                    prep_merge_samples_from_h5ad.out.h5ad.collect(),
                    params.anndata_compression_opts
            )
            file__anndata_merged = merge_samples_from_h5ad.out.anndata
            file__cells_filtered = merge_samples_from_h5ad.out.cells_filtered

        }else if (mode == 'barcodes'){
            log.info """---Using barcode input---"""
            // need a metrics summary for each of them - for 10x and cellbender can use the metrics as an input
            // for deconvolution the default output -i.e - /lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/mo11_work/outputs/test_qc/inputs/file_metadata.tsv
            // format - experiment_id \t celldata1 \celldata2
            prep_merge_samples(channel__file_paths_10x)
            file_metadata.view()
            merge_samples(
                params.output_dir,
                params.input_data_table,
                file_metadata,
                params.file_sample_qc,
                params.file_cellmetadata,
                params.metadata_key_column.value,
                prep_merge_samples.out.barcodes.collect(),
                prep_merge_samples.out.features.collect(),
                prep_merge_samples.out.matrix.collect(),
                params.anndata_compression_opts
            )
            file__anndata_merged = merge_samples.out.anndata
            file__cells_filtered = merge_samples.out.cells_filtered
        }
    emit:
        file__anndata_merged
        file__cells_filtered

}
