include {SCRUBLET} from '../scrublet/main'
include {make_cellmetadata_pipeline_input} from './functions.nf'

workflow MULTIPLET {
    take:
        channel__file_paths_10x
        expected_multiplet_rate
        n_simulated_multiplet
        multiplet_threshold_method
        scale_log10
    
    main:
        // Identify multiplets using scrublet.

        log.info("expected_multiplet_rate: ${params.sample_qc.cell_filters.filter_multiplets.expected_multiplet_rate}")
        // log.info("output_dir: ${output_dir}")
        log.info("n_simulated_multiplet: ${params.sample_qc.cell_filters.filter_multiplets.n_simulated_multiplet}")
        log.info("multiplet_threshold_method: ${params.sample_qc.cell_filters.filter_multiplets.multiplet_threshold_method}")
        log.info("scale_log10: ${params.sample_qc.cell_filters.filter_multiplets.scale_log10}")

        // channel__file_paths_10x.subscribe { println "channel__file_paths_10x: $it" }

        SCRUBLET(
            // output_dir,
            channel__file_paths_10x,
            expected_multiplet_rate,
            n_simulated_multiplet,
            multiplet_threshold_method,
            scale_log10
        )
        
        // Generate input file for merge based in multiplets
        make_cellmetadata_pipeline_input(
            // output_dir,
            SCRUBLET.out.multiplet_calls_published.collect()
        )
    emit:
        // Return merged input data file.
        // outdir = make_cellmetadata_pipeline_input.out.outdir
        file__cellmetadata = make_cellmetadata_pipeline_input.out.file__cellmetadata
        multiplet_calls = SCRUBLET.out.multiplet_calls 
        scrublet_paths = SCRUBLET.out.scrublet_paths
}

