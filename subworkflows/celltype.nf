
include {AZIMUTH;REMAP_AZIMUTH} from "$projectDir/modules/nf-core/modules/azimuth/main"
include {CELLTYPIST} from "$projectDir/modules/nf-core/modules/celltypist/main"
include {SPLIT_BATCH_H5AD} from "$projectDir/modules/nf-core/modules/split_batch_h5ad/main"
include {KERAS_CELLTYPE} from "$projectDir/modules/nf-core/modules/keras_celltype/main"
include {CELLTYPE_FILE_MERGE} from "$projectDir/modules/nf-core/modules/cell_type_assignment/functions"
include {SCPRED} from "$projectDir/modules/nf-core/modules/scpred/main"
include {  DSB } from '../modules/nf-core/modules/citeseq/main'
include { CONVERT_MTX_TO_H5AD } from "$projectDir/modules/local/convert_h5ad_to_mtx/main"

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
        
        // SPLIT_BATCH_H5AD.out.sample_file
        //     .splitCsv(header: true, sep: "\t", by: 1)
        //     .map{row -> tuple(row.experiment_id, file(row.h5ad_filepath))}.set{ch_experiment_filth5}

        // SPLIT_BATCH_H5AD.out.az_sample_file
        //     .splitCsv(header: true, sep: "\t", by: 1)
        //     .map{row -> tuple(row.experiment_id, file(row.h5ad_filepath))}.set{az_ch_experiment_filth5}


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
        
        CELLTYPE_FILE_MERGE(az_out.collect().unique(),
                            ct_out.collect().unique(),
                            all_extra_fields2.collect().unique(),
                            az_ch_experiment_filth5) 

        celltype_assignments=CELLTYPE_FILE_MERGE.out.celltype_assignments

    emit:
        celltype_assignments
}