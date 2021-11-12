include { generate_final_UMAPS } from '../modules/nf-core/modules/umap/functions'

workflow data_handover{
    take:
        predicted_celltypes
        outdir
        file__anndata_merged
        file__cellranger_raw_files_table_tsv
        file__cellbender_files_table_tsv
        file__deconv_files_table_tsv
        multiplet_calls
        deconvolution_path
        qc_output_dir
        outdir_ad
        anndata
    main:
        log.info 'running data handover'

        generate_final_UMAPS(gather_handover_dataset.out.adata_celltypes,anndata,outdir)
}