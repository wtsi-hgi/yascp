
// Load base.config by default for all pipelines - typically included in the nextflow config.
include { prepare_inputs } from './prepare_inputs.nf'
include { main_deconvolution } from './main_deconvolution.nf'
include { CELLTYPIST } from '../modules/nf-core/modules/celltypist/main'

workflow deconvolution {
    // prepare input channels, depending on which input mode was chosen:
    //if (! file(params.input_data_table).isEmpty()) {
	prepare_inputs(Channel.fromPath(params.input_data_table, followLinks: true, checkIfExists: true))

    // run main deconvolution pipeline on prepared input channels:
    main_deconvolution(prepare_inputs.out.ch_experiment_bam_bai_barcodes,
		       prepare_inputs.out.ch_experiment_npooled,
		       prepare_inputs.out.ch_experiment_filth5,
		       prepare_inputs.out.ch_experiment_donorsvcf_donorslist)

	if (params.celltypist.run) {
		// read filtered barcodes straight from cellranger outputs
		log.info "params.celltypist.path_filtered_h5: ${params.celltypist.path_filtered_h5}"

		channel.fromPath(params.celltypist.path_filtered_h5)
				.splitCsv(header: true, sep: params.input_tables_column_delimiter)
			.map{row->tuple(row.experiment_id, row.data_path_10x_format.replaceFirst(/${params.replace_in_path_from}/, params.replace_in_path_to))}
			.set{ch_experiment_filth5}
		// to finish implementing
		channel.fromList(params.celltypist.models)
			.set{ch_celltypist_models}
		// ch_experiment_filth5.combine(ch_celltypist_models).view()
		CELLTYPIST(ch_experiment_filth5.combine(ch_celltypist_models))
	}

}
