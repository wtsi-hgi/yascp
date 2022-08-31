process CREATE_ARTIFICIAL_BAM_CHANNEL_process{
    
    // This will be a csv file that will be later converted to the channel as per data emmited from main_deconvolution channel
    scratch false      // use tmp directory
    label 'process_low'
    input:
        path(outdir)
        path(fech_folder)

    output:
        val(outdir, emit: ch_experiment_bam_bai_barcodes)

    script:
        """
            echo 'lets do this'
        """
        
}



workflow CREATE_ARTIFICIAL_BAM_CHANNEL{
    take:
        input_channel
    main:
        // CREATE_ARTIFICIAL_BAM_CHANNEL_process(params.outdir,input_channel)
        // ch_experiment_bam_bai_barcodes=CREATE_ARTIFICIAL_BAM_CHANNEL_process.out.ch_experiment_bam_bai_barcodes

        input_channel
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.experiment_id, "${row.data_path_10x_format}/possorted_genome_bam.bam","${params.outdir}/deconvolution/vireo/${row.experiment_id}/donor_ids.tsv")}
            .set{ch_experiment_bam_bai_barcodes}
        

        input_channel
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.experiment_id, "${params.outdir}/gtmatch/${row.experiment_id}/${row.experiment_id}_gt_donor_assignments.csv")}
            .set{ch_poolid_csv_donor_assignments}            

    emit:
        ch_experiment_bam_bai_barcodes
        ch_poolid_csv_donor_assignments

}