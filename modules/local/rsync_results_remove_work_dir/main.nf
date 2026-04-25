
process RSYNC_RESULTS_REMOVE_WORK_DIR {
    // This process rsyncs the results to a directory to extract the data from symlings before 
    // deleting work dir in the onComplete statement of pipeline
    
    tag "${pool_id}"
    label 'process_low'
    publishDir  path: "${params.outdir}/deconvolution/concordances/${pool_id}",
                mode: "copy",
                overwrite: "true"

    input:
        path(outdir)
        path(workdir)

    script:

        """
        dir=\$(readlink -f $outdir)
        work_dir=\$(readlink -f $workdir)
        rsync -vrL --update \$dir/ \${dir}_rsync
        echo \$work_dir
        echo \${dir}_rsync
        #rm -r \$work_dir
        rm -r \${dir}
        """
}