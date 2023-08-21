
process RSYNC_RESULTS_REMOVE_WORK_DIR {

    tag "${pool_id}"
    label 'process_low'
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
    //     //// container "https://yascp.cog.sanger.ac.uk/public/singularity_images/mercury_scrna_deconvolution_latest.img"
    // } else {
    //     container "mercury/scrna_deconvolution:62bd56a"
    // }

    publishDir  path: "${params.outdir}/concordances/${pool_id}",
                mode: "copy",
                overwrite: "true"

    input:
        path(outdir)
        path(workdir)

    script:

        """
        dir=\$(readlink -f $outdir)
        work_dir=\$(readlink -f $workdir)
        rsync -vrL \$dir \${dir}_rsync
        echo \$work_dir
        echo \${dir}_rsync
        rm -r \$work_dir
        rm -r \${dir}
        """
}