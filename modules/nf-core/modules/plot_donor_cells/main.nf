
process PLOT_DONOR_CELLS {

    tag "${sample_donor_summary_tsv}"
    
    label 'process_low'
    
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/mercury_scrna_deconvolution_62bd56a-2021-12-15-4d1ec9312485.sif"
        //// container "/software/hgi/containers/mercury_scrna_deconvolution_latest.img"
    } else {
        container "mercury/scrna_deconvolution:62bd56a"
    }

    publishDir "${params.outdir}/plots/", mode: "${params.plot_donor_ncells.copy_mode}", overwrite: true,
	  saveAs: {filename -> filename.indexOf(".pdf") > 0 ? filename.replaceFirst("outputs/","") : "$filename"}
    
    when: 
    params.plot_donor_ncells.run

    input: 
    path(sample_donor_summary_tsv)

    output: 
    path("outputs/*.pdf"), emit: sample_pdf

    script:
    """
        python plot_donor_ncells.py \\
        --output_dir \$PWD/outputs \\
        --sample_donor_summary_tsv ${sample_donor_summary_tsv} \\
        --plotnine_dpi ${params.plot_donor_ncells.plotnine_dpi}
    """
}
