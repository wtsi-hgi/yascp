
// Load base.config by default for all pipelines - typically included in the nextflow config.
include { CELLBENDER } from '../modules/nf-core/modules/cellbender/main'
// Modules to include.
include {
    wf__multiplet;
} from "../modules/nf-core/modules/multiplet/main"


workflow qc {
    main:
        log.info "running QC metrics"
        // log.info params.metadata_key_column.value
        if (run_multiplet_filters) {
	      log.info "Running multiplet filters."
        } else {
            file_cellmetadata = file(params.file_cellmetadata)
            multiplet_calls = null
        }

        // The mode of input may change - Conventional and Subclustering

}
