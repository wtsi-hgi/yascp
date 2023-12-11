process SCPRED{
    tag "${samplename}"    
    label 'process_medium'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "/software/hgi/containers/yascp/cell_classification.sif "
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

    output:
        path('*.RDS')
        // path(celltype_table, emit:predicted_celltypes)
        // path(celltype_table, emit:predicted_celltype_labels)
        // path "ncells_by_type_barplot.pdf"
        // path "query_umap.pdf"
        // path "prediction_score_umap.pdf"
        // path "prediction_score_vln.pdf"
        // path "mapping_score_umap.pdf"
        // path "mapping_score_vln.pdf"

    script:
    
    outdir = "${outdir_prev}/azimuth"
    // output file prefix: strip random hex number form beginning of file name
    outfil_prfx = "${file_h5ad_batch}".minus(".h5ad")
    //outfil_prfx = "${file_h5ad_batch}".minus(".h5ad")
    celltype_table = "${outfil_prfx}_predicted_celltype_l2.tsv.gz"
    """
        scpred_map_hiers.R --file ./${file_h5ad_batch}
    """
}
