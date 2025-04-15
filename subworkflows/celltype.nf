
include {AZIMUTH;REMAP_AZIMUTH} from "$projectDir/modules/nf-core/modules/azimuth/main"
include {CELLTYPIST} from "$projectDir/modules/nf-core/modules/celltypist/main"
include {SPLIT_BATCH_H5AD} from "$projectDir/modules/nf-core/modules/split_batch_h5ad/main"
include {KERAS_CELLTYPE} from "$projectDir/modules/nf-core/modules/keras_celltype/main"
// include {CELLTYPE_FILE_MERGE} from "$projectDir/modules/nf-core/modules/cell_type_assignment/functions"
include {SCPRED} from "$projectDir/modules/nf-core/modules/scpred/main"
include {  DSB } from '../modules/nf-core/modules/citeseq/main'
include { CONVERT_MTX_TO_H5AD } from "$projectDir/modules/local/convert_h5ad_to_mtx/main"

process CELLTYPE_FILE_MERGE{
    tag "${samplename}"    
    label 'process_high'
    publishDir  path: "${params.outdir}/celltype_assignemt/",
            saveAs: {filename ->
                    if (filename.contains("adata.h5ad")) {
                        null
                    } else {
                        filename
                    }
                },
            mode: "${params.copy_mode}",
            overwrite: "true"  

    publishDir  path: "${params.outdir}/handover/merged_h5ad/",
            saveAs: {filename ->
                    if (filename.contains("adata.h5ad")) {
                        filename = "2.celltype_anotated_merged.h5ad"
                    } else {
                        null
                    }
                },
            mode: "${params.copy_mode}",
            overwrite: "true"  


    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"
    } else {
       container "${params.yascp_container_docker}"
    }
    output:
        // path('adata.h5ad', emit:file__anndata_merged2)
        path("All_Celltype_Assignments.tsv",emit:celltype_assignments)
        path "tranche_celltype_report.tsv"
        path "donor_celltype_report.tsv"

    input:
        path(azimuth_files)
        path(celltypist_paths)
        path(all_other_paths)
    script:
        def merged_files_outpath = workflow.workDir.toString()
        file(merged_files_outpath).mkdirs()
        def azimuth_files_path = "${merged_files_outpath}/azimuth_files.tsv"
        def celltypist_files_path = "${merged_files_outpath}/celltypist_files.tsv"
        def all_other_files_path = "${merged_files_outpath}/other_files.tsv"

        new File(azimuth_files_path).text = azimuth_files.join("\n")
        new File(celltypist_files_path).text = celltypist_paths.join("\n")

        if ("${all_other_paths}" != 'fake_file.fq') {
            new File(all_other_files_path).text = all_other_paths.join("\n")
            other_paths = "--all_other_paths ${all_other_files_path}"
        } else {
            other_paths = ""
        }

        """
        generate_combined_celltype_anotation_file.py --all_azimuth_files ${azimuth_files_path} --all_celltypist_files ${celltypist_files_path} ${other_paths}
        """

}


workflow celltype{
    
    take:
        file__anndata_merged
        hastag_labels
        
    main:

        log.info '---Splitting the assignment for each batch---'
        // file__anndata_merged.subscribe { println "file__anndata_merged: $it" }


        // Here we may want to not split it and just pass in an entire h5ad file for annotations.
        // We need a combined h5ad file with all donors to perform further data integrations
        file__anndata_merged_post = CONVERT_MTX_TO_H5AD(file__anndata_merged).gex_h5ad
        // SPLIT_BATCH_H5AD(file__anndata_merged_post,params.split_ad_per_bach)
        ch_experiment_filth5 = file__anndata_merged_post 
        az_ch_experiment_filth5 = file__anndata_merged_post


        // Keras celltype assignemt
        if (params.celltype_assignment.run_keras){
            KERAS_CELLTYPE(ch_experiment_filth5,params.celltype_prediction.keras.keras_model,params.celltype_prediction.keras.keras_weights_df) 
            all_extra_fields3 = KERAS_CELLTYPE.out.predicted_celltype_labels.collect()
            all_extra_fields = all_extra_fields3.ifEmpty(Channel.from("$projectDir/assets/fake_file.fq"))
        }else{
            all_extra_fields = Channel.from("$projectDir/assets/fake_file.fq")
        }
        
        // AZIMUTH
        if (params.celltype_assignment.run_azimuth){
            AZIMUTH(az_ch_experiment_filth5,Channel.fromList( params.azimuth.celltype_refsets))
            az_out = AZIMUTH.out.predicted_celltype_labels.collect()
            // REMAP_AZIMUTH(AZIMUTH.out.celltype_tables_all,params.mapping_file)
            // az_out = REMAP_AZIMUTH.out.predicted_celltype_labels.collect()
        }else{
            az_out = Channel.from("$projectDir/assets/fake_file1.fq")
            az_out = az_out.ifEmpty(Channel.from("$projectDir/assets/fake_file1.fq"))
        }
        
        // CELLTYPIST
        if (params.celltype_assignment.run_celltypist){
            Channel.fromList(params.celltypist.models)
                .set{ch_celltypist_models}
            CELLTYPIST(az_ch_experiment_filth5.combine(ch_celltypist_models))
            ct_out2 = CELLTYPIST.out.predicted_labels.collect()
            ct_out = ct_out2.ifEmpty(Channel.from("$projectDir/assets/fake_file2.fq"))
        }else{
            ct_out = Channel.from("$projectDir/assets/fake_file2.fq")
        }

        // // SCPRED
        if (params.celltype_assignment.run_scpred){
            SCPRED(ch_experiment_filth5,params.scpred.reference)
            sc_out2 = SCPRED.out.predicted_celltype_labels.collect()
            sc_out = sc_out2.ifEmpty(Channel.of())
        }else{
            sc_out = Channel.of()
        }        
        all_extra_fields2 = all_extra_fields.mix(sc_out).mix(hastag_labels)
        CELLTYPE_FILE_MERGE(az_out.collect().unique(),ct_out.collect().unique(),all_extra_fields2.collect().unique()) 
        
        celltype_assignments=CELLTYPE_FILE_MERGE.out.celltype_assignments

    emit:
        celltype_assignments
}