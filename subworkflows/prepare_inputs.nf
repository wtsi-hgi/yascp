nextflow.enable.dsl=2
include { from_barcodes } from './prepare_inputs/from_barcodes.nf'

workflow prepare_inputs {
	// this workflow processes the outputs from cellbender to perform the data preparation
	// TODO: Curently it has been left as close as possible to the original imput of the celbender pipeline, however the sub workflow in  the prepeare_inputs folder is redundant and needs to be cleaned up. 
	// TODO: We shuld also add an input check to make sure that correct format is provided, as per template in locals/inputs provided by nfCore.

    take: channel_input_data_table
    main:

		log.info " ----Prepearing the inputs ---"
		from_barcodes(channel_input_data_table)
		ch_experiment_bam_bai_barcodes = from_barcodes.out.ch_experiment_bam_bai_barcodes
		ch_experiment_npooled = from_barcodes.out.ch_experiment_npooled
		ch_experiment_filth5 = from_barcodes.out.ch_experiment_filth5
		ch_experiment_donorsvcf_donorslist = from_barcodes.out.ch_experiment_donorsvcf_donorslist
		ch_experimentid_paths10x_raw = from_barcodes.out.ch_experimentid_paths10x_raw
		ch_experimentid_paths10x_filtered=from_barcodes.out.ch_experimentid_paths10x_filtered
		channel__file_paths_10x=from_barcodes.out.channel__file_paths_10x
		channel__metadata=from_barcodes.out.channel__metadata

    emit:
		ch_experiment_bam_bai_barcodes
		ch_experiment_npooled
		ch_experiment_filth5
		ch_experiment_donorsvcf_donorslist
        ch_experimentid_paths10x_raw
        ch_experimentid_paths10x_filtered
        channel__file_paths_10x
        channel__metadata
}
