
// Load base.config by default for all pipelines - typically included in the nextflow config.
include { CELLBENDER } from '../modules/nf-core/modules/cellbender/main'

workflow cellbender {
    take:
        ch_experimentid_paths10x_raw
		ch_experimentid_paths10x_filtered
        channel__metadata
    main:
        log.info params.input_data_table
        log.info """---Running Cellbender pipeline ---"""
        CELLBENDER(ch_experimentid_paths10x_raw,ch_experimentid_paths10x_filtered,channel__metadata)
        // results_list = CELLBENDER.out.results_list
        cellbender_path=CELLBENDER.out.cellbender_path
    emit:
        // results_list
        cellbender_path
}
