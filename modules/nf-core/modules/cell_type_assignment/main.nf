
include {AZIMUTH} from '../azimuth/main'
include {CELLTYPIST} from '../celltypist/main'
include {SPLIT_BATCH_H5AD} from '../split_batch_h5ad/main'
include {CELLTYPE_FILE_MERGE} from './functions'

workflow CELL_TYPE_ASSIGNEMT{
    
    take:
        file__anndata_merged
        file__cells_filtered
        
    main:
        SPLIT_BATCH_H5AD(file__anndata_merged)
        SPLIT_BATCH_H5AD.out.sample_file.view()
        SPLIT_BATCH_H5AD.out.sample_file
            .splitCsv(header: true, sep: "\t", by: 1)
            .map{row -> tuple(row.experiment_id, file(row.h5ad_filepath))}.set{ch_experiment_filth5}
        Channel.fromList(params.celltypist.models)
            .set{ch_celltypist_models}
        ch_experiment_filth5.combine(ch_celltypist_models)

        SPLIT_BATCH_H5AD.out.files_anndata_batch.flatMap().set{ch_batch_files}
        AZIMUTH(params.output_dir,ch_batch_files)
        CELLTYPIST(ch_experiment_filth5.combine(ch_celltypist_models))
        // Add all the assignments to the file__anndata_merged and return it for downstraeam analysis
        CELLTYPE_FILE_MERGE(AZIMUTH.out.predicted_celltype_labels.collect(),CELLTYPIST.out.sample_predicted_labels_csv.collect(),file__anndata_merged)
        file__anndata_merged2=CELLTYPE_FILE_MERGE.out.file__anndata_merged2

    emit:
        file__anndata_merged2
}