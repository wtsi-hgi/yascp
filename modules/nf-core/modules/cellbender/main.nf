// Modules to include.
include {
    cellbender__rb__get_input_cells;
    cellbender__remove_background;
    cellbender__remove_background__qc_plots;
    cellbender__remove_background__qc_plots_2;
    cellbender__gather_qc_input;cellbender__preprocess_output;
} from "./functions.nf"

// Set default parameters.
outdir           = "${params.output_dir}/nf-preprocessing"
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
        channel__metadata
        
    main:
  
        ch_experimentid_paths10x_raw.map{row -> tuple(
            row[0],
            file("${row[1]}/barcodes.tsv.gz"),
            file("${row[1]}/features.tsv.gz"),
            file("${row[1]}/matrix.mtx.gz")
        )}.set{channel__file_paths_10x}
    
        outdir =  outdir+'/cellbender'
        
        // here pass in the number of cells detected by cellranger/ 
        if (params.cellbender_rb.estimate_params_umis.value.method_estimate_ncells=='expected'){
            channel__metadata.splitCsv(header: true, sep: "\t", by: 1).map{row -> tuple(
                row.experiment_id,
                row.Estimated_Number_of_Cells,
            )}.set{ncells_cellranger_pre}
        }else{
            channel__metadata.splitCsv(header: true, sep: "\t", by: 1).map{row -> tuple(
                row.experiment_id,
                '0',
            )}.set{ncells_cellranger_pre}
        }


        
        ncells_cellranger_pre.join(ch_experimentid_paths10x_raw, remainder: false).set{post_ncells_cellranger} 

        post_ncells_cellranger.map{row -> tuple(row[0], row[1])}.filter{ it[2] == null }.set{ncells_cellranger}
                

        channel__file_paths_10x.combine(ncells_cellranger, by: 0).set{channel__file_paths_10x_with_ncells}
       
        cellbender__rb__get_input_cells(
            outdir,
            channel__file_paths_10x_with_ncells,
            params.cellbender_rb.estimate_params_umis.value,
        )
        
        // Correct counts matrix to remove ambient RNA
        cellbender__rb__get_input_cells.out.cb_input.subscribe { println "cellbender__rb__get_input_cells: $it" }
        cellbender__remove_background(
            outdir,
            cellbender__rb__get_input_cells.out.cb_input,
            params.cellbender_rb.epochs.value,
            params.cellbender_rb.learning_rate.value,
            params.cellbender_rb.zdim.value,
            params.cellbender_rb.zlayers.value,
            params.cellbender_rb.low_count_threshold.value,
            params.cellbender_rb.fpr.value
        )

        cellbender__preprocess_output(
            cellbender__remove_background.out.cleanup_input,
            cellbender__remove_background.out.cb_plot_input,
            cellbender__remove_background.out.experimentid_outdir_cellbenderunfiltered_expectedcells_totaldropletsinclude,
        )

        // Make some basic plots
        cellbender__remove_background__qc_plots(
            cellbender__preprocess_output.out.cb_plot_input
        )

        cellbender__preprocess_output.out.experimentid_outdir_cellbenderunfiltered_expectedcells_totaldropletsinclude
            .combine(ch_experimentid_paths10x_raw, by: 0)
            .combine(ch_experimentid_paths10x_filtered, by: 0)
            .combine(Channel.from("${params.cellbender_rb.fpr.value}"
                    .replaceFirst(/]$/,"")
                    .replaceFirst(/^\[/,"")
                    .split()))
            .set{input_channel_qc_plots_2}
            
            cellbender__remove_background__qc_plots_2(input_channel_qc_plots_2,outdir)
            

        results_list = cellbender__preprocess_output.out.out_paths
        // prepeare the output channel for utilising in the deconvolution instead of barcode input.
        cellbender_path = cellbender__preprocess_output.out.alternative_input

        emit:
            // results_list //results list is not needed to be emited - if done it will wait for all the cellbender stuff to finish.
            cellbender_path

            
}

