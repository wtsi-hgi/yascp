process SCPRED{
    tag "${samplename}"    
    label 'process_medium'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
        container "wtsihgi/nf_scrna_qc_azimuth:d54db9b"
    }

    publishDir  path: "${params.outdir}/celltype_assignemt/scpred/",
            saveAs: {filename -> "${outfil_prfx}_" + filename},
            mode: "${params.copy_mode}",
            overwrite: "true"
    stageInMode 'copy'

    input:
        tuple val(exp_id), path(file_h5ad_batch)
        path(reference)
    when:
        params.celltype_assignment.run_azimuth
    output:
        // path('*.RDS')
        path("${outfil_prfx}___scpred_prediction.tsv"), emit:predicted_celltype_labels

    script:

        outfil_prfx = "${file_h5ad_batch}".minus(".h5ad")
        celltype_table = "${outfil_prfx}_predicted_celltype_l2.tsv.gz"
        """
            scpred_map_hiers.R --file ./${file_h5ad_batch} --reference ${reference}
            ln -s scpred_prediction.tsv ${outfil_prfx}___scpred_prediction.tsv
        """
}
