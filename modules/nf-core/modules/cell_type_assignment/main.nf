
include {AZIMUTH} from '../azimuth/main'
include {CELLTYPIST} from '../celltypist/main'
include {SPLIT_BATCH_H5AD} from '../split_batch_h5ad/main'

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
}