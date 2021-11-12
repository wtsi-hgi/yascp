
def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}


process make_cellmetadata_pipeline_input {
    // Makes a input tsv file for the main pipeline.
    // ------------------------------------------------------------------------
    //tag { output_dir }
    //cache false        // cache results from run

    publishDir  path: "${outdir}",
                saveAs: {filename -> filename.replaceAll("${runid}-", "")},
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        val(outdir_prev)
        path("*multiplet_calls_published.txt")

    output:
        val(outdir, emit: outdir)
        path('file_cellmetadata.tsv', emit: file__cellmetadata)

    script:
        runid = random_hex(16)
        outdir = "${outdir_prev}"
        process_info = "${runid} (runid)"
        process_info = "${process_info}, ${task.cpus} (cpus)"
        process_info = "${process_info}, ${task.memory} (memory)"
        """
        # Note: the default paste delim is tab
        cat *multiplet_calls_published.txt \
            | awk 'BEGIN{print "experiment_id\tdata_path_cellmetadata"}1' \
            > file_cellmetadata.tsv
        """
}