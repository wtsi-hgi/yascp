process AZIMUTH{
    tag "${samplename}"    
    label 'process_medium'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.nf_scrna_qc_azimuth_container}"        
    } else {
        container "wtsihgi/nf_scrna_qc_azimuth:d54db9b"
    }

    publishDir  path: "${params.outdir}/celltype/azimuth/${refset.name}",
            saveAs: {filename -> "${outfil_prfx}_" + filename},
            mode: "${params.copy_mode}",
            overwrite: "true"
    stageInMode 'copy'
    // stageInMode 'copy' is needed because SeuratDisk:::Convert()
    // generates the output file apparently from the absolute path of input file.
    // Symbolic links have the output file written to the link target directory
    // where it cannot be found by the azimuth.R script.

    input:
        val outdir_prev
        tuple val(samplename),path(file_h5ad_batch)
        each refset
        // path(mapping_file)
    output:
        // path(celltype_table, emit:predicted_celltypes)
        tuple(val(outfil_prfx), val(refset.refset), path("*predicted_*.tsv"),emit:celltype_tables_all) 
        path("*predicted_*.tsv"), emit:predicted_celltype_labels
        path "*ncells_by_type_barplot.pdf"
        path "*query_umap.pdf"
        path "*prediction_score_umap.pdf"
        path "*prediction_score_vln.pdf"
        path "*mapping_score_umap.pdf"
        path "*mapping_score_vln.pdf"
        // path("${outfil_prfx}_query.rds"), emit: query_rds

        
    script:
    
    
        // output file prefix: strip random hex number form beginning of file name
        outfil_prfx = "${file_h5ad_batch}".minus(".h5ad")
        //outfil_prfx = "${file_h5ad_batch}".minus(".h5ad")
        if (refset.refset =='PBMC' && params.mapping_file!='' && params.remap_celltypes){

            com="remap_azimuth_l2.R --out_file ${samplename}_predicted_celltype_l2.tsv --mapping ${params.mapping_file} --az_file ${samplename}_predicted_celltype_l2.tsv"
        }else{
            com=""
        }
    
    """ 
        azimuth.R ./${file_h5ad_batch} ${refset.refset} ${refset.annotation_labels} ${samplename}
        echo ${params.mapping_file}
        echo ${params.remap_celltypes}
        ${com}
    """
}

process REMAP_AZIMUTH{
    // This process remaps Azimuth L2 to L1 and L0
    tag "${samplename}"    
    label 'process_low'
   
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.nf_scrna_qc_sif_container}"
    } else {
        container "wtsihgi/nf_scrna_qc_azimuth:d54db9b"
    }

    publishDir  path: "${params.outdir}/celltype/",
            mode: "${params.copy_mode}",
            overwrite: "true"
    stageInMode 'copy'  

    input:
        path(file)
        path(mapping_file)

    output:
        path("remapped__${file}", emit:predicted_celltype_labels)

    script:
        """
            remap_azimuth_l2.py -of remapped__${file} -m ${mapping_file} -az ${file}
        """

}
