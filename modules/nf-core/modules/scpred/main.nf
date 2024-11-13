process SCPRED{
    tag "${samplename}"    
    label 'process_medium'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.cell_classification_container}"
    } else {
        container "wtsihgi/nf_scrna_qc_azimuth:d54db9b"
    }

    publishDir  path: "${params.outdir}/celltype/scpred/",
            saveAs: {filename -> "${outfil_prfx}_" + filename},
            mode: "${params.copy_mode}",
            overwrite: "true"
    stageInMode 'copy'

    input:
        val outdir_prev
        path file_h5ad_batch
    when:
        params.celltype_assignment.run_azimuth
    output:
        path('*.RDS')
        path("${outfil_prfx}__scpred_prediction.tsv"), emit:predicted_celltype_labels

    script:
    
    outdir = "${outdir_prev}/azimuth"
    // output file prefix: strip random hex number form beginning of file name
    outfil_prfx = "${file_h5ad_batch}".minus(".h5ad")
    //outfil_prfx = "${file_h5ad_batch}".minus(".h5ad")
    celltype_table = "${outfil_prfx}_predicted_celltype_l2.tsv.gz"
    """
        scpred_map_hiers.R --file ./${file_h5ad_batch}
        ln -s scpred_prediction.tsv ${outfil_prfx}__scpred_prediction.tsv
    """
}
