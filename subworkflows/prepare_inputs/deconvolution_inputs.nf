nextflow.enable.dsl=2
workflow DECONV_INPUTS{
    // This is a function that prpeares the inputs for deconvolution that is utilised twice.
    take:
        cellbender_path
        prepare_inputs
    main:
        cellbender_path
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.experiment_id, file(row.data_path_10x_format))}
            .set{ch_experiment_filth5} // this channel is used for task 'split_donor_h5ad'
        prepare_inputs.out.ch_experiment_bam_bai_barcodes.map { experiment, bam, bai, barcodes -> tuple(experiment,
                            bam,
                            bai)}.set{pre_ch_experiment_bam_bai_barcodes}
        cellbender_path
            .splitCsv(header: true, sep: params.input_tables_column_delimiter).map{row->tuple(row.experiment_id, file("${row.data_path_10x_format}".replaceFirst("results","${params.outdir}")+'/barcodes.tsv.gz'))}.set{barcodes}
        channel__file_paths_10x= cellbender_path
            .splitCsv(header: true, sep: params.input_tables_column_delimiter).map{row->tuple(row.experiment_id,
                                                    file("${row.data_path_10x_format}".replaceFirst("results","${params.outdir}")+'/barcodes.tsv.gz'),
                                                    file("${row.data_path_10x_format}".replaceFirst("results","${params.outdir}")+'/features.tsv.gz'),
                                                    file("${row.data_path_10x_format}".replaceFirst("results","${params.outdir}")+'/matrix.mtx.gz'))}
        pre_ch_experiment_bam_bai_barcodes.combine(barcodes, by: 0).set{ch_experiment_bam_bai_barcodes}
    emit:
        channel__file_paths_10x
        ch_experiment_bam_bai_barcodes
        ch_experiment_filth5
}