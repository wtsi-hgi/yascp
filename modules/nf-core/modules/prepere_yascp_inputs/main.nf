process YASCP_INPUTS {
    // Takes annData object with PCs and returns plots
    // ------------------------------------------------------------------------
    //cache false        // cache results from run
    scratch false      // use tmp directory

    label 'process_low'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif"
        //// container "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/singularity_images/nf_qc_cluster_2.4.img"
        
    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    input:
        path(input_file)

    // the output for this is a correct format input files as per cb 6.1
    output:
        path("input_file_corectly_formatted.tsv"), emit: input_file_corectly_formatted

    script:
        
        """
        link_files.py --input_table ${input_file} \
        """
}