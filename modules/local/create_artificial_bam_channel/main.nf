
workflow CREATE_ARTIFICIAL_BAM_CHANNEL{
    take:
        input_channel
    main:

        input_channel
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.experiment_id, "${row.data_path_10x_format}/possorted_genome_bam.bam","${params.outdir}/deconvolution/vireo_raw/${row.experiment_id}/donor_ids.tsv")}
            .set{bam_split_channel}

        input_channel
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.experiment_id, "${row.data_path_10x_format}/possorted_genome_bam.bam","${row.data_path_10x_format}/possorted_genome_bam.bam.bai","${params.outdir}/nf-preprocessing/cellbender/${row.experiment_id}/cellbender-epochs_250__learnrt_0pt000005__zdim_100__zlayer_500__lowcount_10/cellbender-FPR_${params.cellbender_resolution_to_use}-filtered_10x_mtx/barcodes.tsv.gz")}
            .set{ch_experiment_bam_bai_barcodes}        

        input_channel
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.experiment_id, "${params.outdir}/deconvolution/gtmatch/${row.experiment_id}/${row.experiment_id}_gt_donor_assignments.csv")}
            .set{ch_poolid_csv_donor_assignments}            

    emit:
        bam_split_channel
        ch_poolid_csv_donor_assignments
        ch_experiment_bam_bai_barcodes

}