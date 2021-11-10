
// Load base.config by default for all pipelines - typically included in the nextflow config.
// Modules to include.
include {MULTIPLET} from "../modules/nf-core/modules/multiplet/main"
include {OUTLIER_FILTER} from "../modules/nf-core/modules/outlier_filter/main"

workflow qc {
    take:
        channel__file_paths_10x
        file__anndata_merged
        file__cells_filtered
    main:
        log.info "--- Running QC metrics --- "
        


        // // DETECTING MULTIPLETS
        // if (params.run_multiplet) {
	    //   log.info "Running multiplet filters."
          
        //   MULTIPLET(
        //     params.output_dir,
        //     channel__file_paths_10x,
        //     params.sample_qc.cell_filters.filter_multiplets.expected_multiplet_rate,
        //     params.sample_qc.cell_filters.filter_multiplets.n_simulated_multiplet,
        //     params.sample_qc.cell_filters.filter_multiplets.multiplet_threshold_method,
        //     params.sample_qc.cell_filters.filter_multiplets.scale_log10
        //     )
        //     file_cellmetadata = MULTIPLET.out.file__cellmetadata
        //     multiplet_calls = MULTIPLET.out.multiplet_calls

        // } else {
        //     file_cellmetadata = file(params.file_cellmetadata)
        //     multiplet_calls = null
        // }


        //FILTERING OUTLIER CELLS
        if (params.sample_qc.cell_filters.filter_outliers.run_process) {
            log.info """---Running automatic outlier cell filtering.----"""
            OUTLIER_FILTER(
                params.output_dir,
                file__anndata_merged,
                file__cells_filtered,
                params.sample_qc.cell_filters.filter_outliers.metadata_columns,
                params.sample_qc.cell_filters.filter_outliers.method,
                params.sample_qc.cell_filters.filter_outliers.outliers_fraction,
                params.sample_qc.cell_filters.filter_outliers.max_samples,
                params.anndata_compression_opts
            )
            file__anndata_merged = OUTLIER_FILTER.out.anndata
            file__cells_filtered = OUTLIER_FILTER.out.cells_filtered
        }


        // The mode of input may change - Conventional and Subclustering
        //  if input is h5ad merge all the h5ads, else convert the 10x to h5ad and merge them.
}
