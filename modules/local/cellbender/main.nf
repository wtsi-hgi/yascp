// Modules to include.
include {
    cellbender__rb__get_input_cells;
    cellbender__remove_background;
    cellbender__remove_background__qc_plots;
    cellbender__remove_background__qc_plots_2;
    cellbender__gather_qc_input;cellbender__preprocess_output;
} from "./functions.nf"

// Set default parameters.
outdir           = "${params.outdir}/preprocessing"

workflow CELLBENDER {
    take:
        ch_experimentid_paths10x_raw
		    ch_experimentid_paths10x_filtered
        channel__metadata
        
    main:
        Channel.empty().set { ch_versions }
        ch_experimentid_paths10x_raw.map{row -> tuple(
            row[0],
            file("${row[1]}/barcodes.tsv.gz"),
            file("${row[1]}/features.tsv.gz"),
            file("${row[1]}/matrix.mtx.gz")
        )}.set{channel__file_paths_10x}

        ch_experimentid_paths10x_raw.map{row -> 
            row[0]}.set{experiment_id_in}

        outdir =  outdir+'/cellbender'

        bc_in = Channel.fromList( params.cellbender_rb.per_sample_thresholds)
        bc_in.map{row -> tuple(
            row.name,
            row.total_droplets_included == "" ? '0' : row.total_droplets_included ,
            row.expected_cells == "" ? '0' : row.expected_cells,
        )}.set{ncells_cellranger_pre}


        ch_experimentid_paths10x_raw.join(ncells_cellranger_pre, by: [0], remainder: true).set{post_ncells_cellranger} 
         

        post_ncells_cellranger
            .map { row ->
                if (row[2] == null) {
                    tuple(row[0], row[1], 0, 0)
                } else {
                    tuple(row[0], row[1], row[2], row[3])
                }
            }
            .set { ncells_cellranger }
        // ncells_cellranger.subscribe { println "ncells_cellranger: $it" }

        channel__file_paths_10x.combine(ncells_cellranger, by: 0).set{channel__file_paths_10x_with_ncells}
       

        // channel__file_paths_10x_with_ncells.subscribe { println "channel__file_paths_10x_with_ncells: $it" }

        cellbender__rb__get_input_cells(
            outdir,
            channel__file_paths_10x_with_ncells,
            params.cellbender_rb.estimate_params_umis.value,
        )
        
        // Correct counts matrix to remove ambient RNA
        // Some samples may fail with the defaults. Hence here we allow for a changes to be applied. 

        epochs_to_use = params.cellbender_rb.epochs.value
        learning_rate_to_use = params.cellbender_rb.learning_rate.value
        zdims_to_use = params.cellbender_rb.zdim.value
        zlayers_to_use = params.cellbender_rb.zlayers.value
        low_count_threshold_to_use = params.cellbender_rb.low_count_threshold.value


        bc_in = Channel.fromList( params.cellbender_rb.per_sample_thresholds)
        bc_in.map{row -> tuple(
            row.name,
            row.low_count_threshold == "" ? params.cellbender_rb.low_count_threshold.value : row.low_count_threshold ,
            row.epochs == "" ? params.cellbender_rb.epochs.value : row.epochs,
            row.learning_rate == "" ? params.cellbender_rb.learning_rate.value : row.learning_rate,
            row.zdim == "" ? params.cellbender_rb.zdim.value : row.zdim,
            row.zlayers == "" ? params.cellbender_rb.zlayers.value : row.zlayers,
        )}.set{channel__f}

        // Now we figure out the ones that are not covered by the individual definitions but are in the input panels needed to process
        // cellbender__rb__get_input_cells.out.cb_input.subscribe { println "cellbender__rb__get_input_cells: $it" }

        cellbender__rb__get_input_cells.out.cb_input.join(channel__f, by: [0], remainder: true).set{post_ch_experimentid_paths10x_filtered}
        post_ch_experimentid_paths10x_filtered.filter{ it[8] == null }.map{row -> tuple(row[0])}.set{not_defined}

        not_defined.map{row -> tuple(
            row[0],
            params.cellbender_rb.low_count_threshold.value,
            params.cellbender_rb.epochs.value,
            params.cellbender_rb.learning_rate.value,
            params.cellbender_rb.zdim.value,
            params.cellbender_rb.zlayers.value,
        )}.set{channel__g}
        
        channel__combo =channel__g.concat(channel__f)

        cellbender__rb__get_input_cells.out.cb_input.join(channel__combo, by: [0], remainder: false).set{cellbender_ambient_rna_input}
        
        filteredChan = cellbender_ambient_rna_input
        .filter { tuple ->
            !params.cellbender_ignore_list.contains(tuple[0])
        }


        cellbender__remove_background(
            outdir,
            filteredChan,
            params.cellbender_rb.fpr.value
        )
        ch_versions = ch_versions.mix(cellbender__remove_background.out.versions)

        cellbender__preprocess_output(
            cellbender__remove_background.out.cleanup_input,
            cellbender__remove_background.out.cb_plot_input,
            cellbender__remove_background.out.experimentid_outdir_cellbenderunfiltered_expectedcells_totaldropletsinclude,
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
        cellbender_path_raw = cellbender__preprocess_output.out.alternative_input_raw
        cellbender_downstream = cellbender__remove_background.out.cb_to_use_downstream
        emit:
            // results_list //results list is not needed to be emited - if done it will wait for all the cellbender stuff to finish.
            cellbender_path
            cellbender_downstream
            cellbender_path_raw
            cellbender_versions = ch_versions

            
}

