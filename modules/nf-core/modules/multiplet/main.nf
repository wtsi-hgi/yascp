include {SCRUBLET} from '../scrublet/main'
include {make_cellmetadata_pipeline_input} from './functions.nf'

workflow MULTIPLET {
    take:
        output_dir
        channel__file_paths_10x
        expected_multiplet_rate
        n_simulated_multiplet
        multiplet_threshold_method
        scale_log10
    
    main:
        // Identify multiplets using scrublet.
        SCRUBLET(
            output_dir,
            channel__file_paths_10x,
            expected_multiplet_rate,
            n_simulated_multiplet,
            multiplet_threshold_method,
            scale_log10
        )
        
        // Generate input file for merge based in multiplets
        make_cellmetadata_pipeline_input(
            output_dir,
            SCRUBLET.out.multiplet_calls_published.collect()
        )
    emit:
        // Return merged input data file.
        outdir = make_cellmetadata_pipeline_input.out.outdir
        file__cellmetadata = make_cellmetadata_pipeline_input.out.file__cellmetadata
        multiplet_calls = SCRUBLET.out.multiplet_calls 
}

