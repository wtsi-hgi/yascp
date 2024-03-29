include { SCRUBLET } from "$projectDir/modules/nf-core/modules/scrublet/main"
include { doubletdetection } from "$projectDir/modules/nf-core/modules/doubletdetection/main"

def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

process make_cellmetadata_pipeline_input {
    // Makes a input tsv file for the main pipeline.
    // ------------------------------------------------------------------------
    //cache false        // cache results from run

    publishDir  path: "${params.outdir}/multiplet.method=scrublet",
                saveAs: {filename -> filename.replaceAll("-", "")},
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        path("*multiplet_calls_published.txt")

    output:
        path('file_cellmetadata.tsv', emit: file__cellmetadata)

    script:

        """
        # Note: the default paste delim is tab
        cat *multiplet_calls_published.txt \
            | awk 'BEGIN{print "experiment_id\tdata_path_cellmetadata"}1' \
            > file_cellmetadata.tsv
        """
}

workflow MULTIPLET {
    take:
        channel__file_paths_10x
    
    main:
        // Identify multiplets using scrublet.

        log.info("expected_multiplet_rate: ${params.sample_qc.cell_filters.filter_multiplets.expected_multiplet_rate}")
        log.info("n_simulated_multiplet: ${params.sample_qc.cell_filters.filter_multiplets.n_simulated_multiplet}")
        log.info("multiplet_threshold_method: ${params.sample_qc.cell_filters.filter_multiplets.multiplet_threshold_method}")
        log.info("scale_log10: ${params.sample_qc.cell_filters.filter_multiplets.scale_log10}")

        SCRUBLET(
            channel__file_paths_10x,
            params.sample_qc.cell_filters.filter_multiplets.expected_multiplet_rate,
            params.sample_qc.cell_filters.filter_multiplets.n_simulated_multiplet,
            params.sample_qc.cell_filters.filter_multiplets.multiplet_threshold_method,
            params.sample_qc.cell_filters.filter_multiplets.scale_log10
        )
        
        // doubletdetection(channel__file_paths_10x)

        // Generate input file for merge based in multiplets
        make_cellmetadata_pipeline_input(
            SCRUBLET.out.multiplet_calls_published.collect()
        )

    emit:
        // Return merged input data file.
        // outdir = make_cellmetadata_pipeline_input.out.outdir
        file__cellmetadata = make_cellmetadata_pipeline_input.out.file__cellmetadata
        multiplet_calls = SCRUBLET.out.multiplet_calls 
        scrublet_paths = SCRUBLET.out.scrublet_paths
}

