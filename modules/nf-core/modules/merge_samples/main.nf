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
        mode
    main:
        log.info """---Merging samples in a single h5ad file---"""

        if (mode == 'h5ad'){
            log.info """---Using h5ad input file outputs from deconvolution---"""
            
            channel__file_paths_10x.splitCsv(header: true, sep: "\t", by: 1)
                .map{row -> tuple(row.experiment_id, file(row.h5ad_filepath))}.set{channel__file_paths_10x_paths}
        
            prep_merge_samples_from_h5ad(channel__file_paths_10x_paths)
            
            // This is currently not finished -- currently hard filters happen in the merging part next. We may want to keep this seperate so we can filter the stuff out after celltype assignments.
            // CELL_HARD_FILTERS(prep_merge_samples_from_h5ad.out.h5ad,params.hard_filters_file) //ad with all cells goes in and adata with removed cells comes out.
            merge_samples_from_h5ad(
                    params.outdir,
                    channel__file_paths_10x,
                    file_metadata,
                    params.file_cellmetadata,
                    params.metadata_key_column.value,
                    prep_merge_samples_from_h5ad.out.h5ad.collect(),
                    params.anndata_compression_opts
            )
            file__anndata_merged = merge_samples_from_h5ad.out.anndata
            // file__cells_filtered = merge_samples_from_h5ad.out.cells_filtered

        }else if (mode == 'barcodes'){
            log.info """---Using barcode input---"""
            // need a metrics summary for each of them - for 10x and cellbender can use the metrics as an input
            // for deconvolution the default output -i.e - /lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/mo11_work/outputs/test_qc/inputs/file_metadata.tsv
            // format - experiment_id \t celldata1 \celldata2
            prep_merge_samples(channel__file_paths_10x)
            
            merge_samples(
                params.outdir,
                params.input_data_table,
                file_metadata,
                params.file_cellmetadata,
                params.metadata_key_column.value,
                prep_merge_samples.out.barcodes.collect(),
                prep_merge_samples.out.features.collect(),
                prep_merge_samples.out.matrix.collect(),
                params.anndata_compression_opts
            )
            file__anndata_merged = merge_samples.out.anndata
            // file__cells_filtered = merge_samples.out.cells_filtered
        }
    emit:
        file__anndata_merged
        // file__cells_filtered

}
