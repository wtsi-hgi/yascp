
// Load base.config by default for all pipelines - typically included in the nextflow config.
include { CELLBENDER } from '../modules/nf-core/modules/cellbender/main'

workflow cellbender {
    log.info params.input_data_table
    log.info """--- Running Cellbender pipeline ---"""
    CELLBENDER()
}
