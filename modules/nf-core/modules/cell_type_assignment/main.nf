
include {AZIMUTH} from "$projectDir/modules/nf-core/modules/azimuth/main"
include {CELLTYPIST} from "$projectDir/modules/nf-core/modules/celltypist/main"
include {SPLIT_BATCH_H5AD} from "$projectDir/modules/nf-core/modules/split_batch_h5ad/main"
include {KERAS_CELLTYPE} from "$projectDir/modules/nf-core/modules/keras_celltype/main"
        
include {CELLTYPE_FILE_MERGE} from './functions'

workflow CELL_TYPE_ASSIGNEMT{
    
    take:
        file__anndata_merged
        file__cells_filtered
        
    main:
        // if (params.split_ad_per_bach){
        log.info '---Splitting the assignment for each batch---'

        file__anndata_merged.map{val1 -> tuple('full_ct', val1)}.set{out1}
        // KERAS_CELLTYPE(out1)   

        SPLIT_BATCH_H5AD(file__anndata_merged,params.split_ad_per_bach)
        // SPLIT_BATCH_H5AD.out.sample_file.view()
        // Here we may want to not split it and just pass in an entire h5ad file for annotations.
        SPLIT_BATCH_H5AD.out.sample_file
            .splitCsv(header: true, sep: "\t", by: 1)
            .map{row -> tuple(row.experiment_id, file(row.h5ad_filepath))}.set{ch_experiment_filth5}
        SPLIT_BATCH_H5AD.out.az_sample_file
            .splitCsv(header: true, sep: "\t", by: 1)
            .map{row -> tuple(row.experiment_id, file(row.h5ad_filepath))}.set{az_ch_experiment_filth5}

        SPLIT_BATCH_H5AD.out.files_anndata_batch.flatMap().set{ch_batch_files}

        KERAS_CELLTYPE(ch_experiment_filth5) 
        Channel.fromList(params.celltypist.models)
            .set{ch_celltypist_models}

        ch_experiment_filth5.view()    
        AZIMUTH(params.output_dir,ch_batch_files)
        CELLTYPIST(az_ch_experiment_filth5.combine(ch_celltypist_models))
        
        CELLTYPE_FILE_MERGE(AZIMUTH.out.predicted_celltype_labels.collect(),CELLTYPIST.out.sample_predicted_labels_csv.collect(),KERAS_CELLTYPE.out.predicted_celltype_labels.collect(),file__anndata_merged)
        file__anndata_merged2=CELLTYPE_FILE_MERGE.out.file__anndata_merged2

    emit:
        file__anndata_merged2
}