// Modules to include.
include {
    cellbender__rb__get_input_cells;
    cellbender__remove_background;
    cellbender__remove_background__qc_plots;
    cellbender__remove_background__qc_plots_2;
    cellbender__gather_qc_input;
} from "./core.nf"

// Set default parameters.
params.output_dir           = "nf-preprocessing"
params.help                 = false

params.cellbender_rb = [
    estimate_params_umis: [value: [
        expected_nemptydroplets_umi_cutoff: 0,
        method_estimate_ncells: 'dropletutils::barcoderanks::inflection',
        lower_bound_umis_estimate_ncells: 1000,
        method_estimate_nemptydroplets: 'dropletutils::barcoderanks::inflection',
        lower_bound_umis_estimate_nemptydroplets: 10,
        upper_bound_umis_estimate_nemptydroplets: 100,
        estimate_nemptydroplets_umi_subtract_factor: 25
    ]],
    epochs: [value: [200]],
    learning_rate: [value: [0.001, 0.0001]],
    fpr: [value: [0.01, 0.05]],
]

workflow CELLBENDER {
    take:
        ch_experimentid_paths10x_raw
		ch_experimentid_paths10x_filtered
        
    main:
        channel__file_paths_10x = Channel
        .fromPath(params.input_data_table)
        .splitCsv(header: true, sep: "\t", by: 1)
        .map{row -> tuple(
            row.experiment_id,
            file("${row.data_path_10x_format}/raw_feature_bc_matrix/barcodes.tsv.gz"),
            file("${row.data_path_10x_format}/raw_feature_bc_matrix/features.tsv.gz"),
            file("${row.data_path_10x_format}/raw_feature_bc_matrix/matrix.mtx.gz")
        )}
    

        
        cellbender__rb__get_input_cells(
            params.output_dir,
            channel__file_paths_10x,
            params.cellbender_rb.estimate_params_umis.value
        )
        
        // Correct counts matrix to remove ambient RNA
        cellbender__remove_background(
            params.output_dir,
            cellbender__rb__get_input_cells.out.cb_input,
            params.cellbender_rb.epochs.value,
            params.cellbender_rb.learning_rate.value,
            params.cellbender_rb.zdim.value,
            params.cellbender_rb.zlayers.value,
            params.cellbender_rb.low_count_threshold.value,
            params.cellbender_rb.fpr.value
        )


        // Make some basic plots
        cellbender__remove_background__qc_plots(
            cellbender__remove_background.out.cb_plot_input
        )

        // if (params.cellbender__remove_background__qc_plots_2.run_task) {
 
            log.info "will run 'cellbender__remove_background__qc_plots_2' Nextflow task"
            
            cellbender__remove_background.out.experimentid_outdir_cellbenderunfiltered_expectedcells_totaldropletsinclude
                .combine(ch_experimentid_paths10x_raw, by: 0)
                .combine(ch_experimentid_paths10x_filtered, by: 0)
                .combine(Channel.from("${params.cellbender_rb.fpr.value}"
                        .replaceFirst(/]$/,"")
                        .replaceFirst(/^\[/,"")
                        .split()))
                .set{input_channel_qc_plots_2}
            
            cellbender__remove_background__qc_plots_2(input_channel_qc_plots_2)
            
        // }
    
        cellbender__gather_qc_input(
            params.output_dir,
            cellbender__remove_background.out.results_list.collect()
        )

        results_list = cellbender__remove_background.out.results_list
        results_list.view()

        // prepeare the output channel for utilising in the deconvolution instead of barcode input.
        


        emit:
            results_list
            // TODO here need to emmit the files from cellbender.
}

