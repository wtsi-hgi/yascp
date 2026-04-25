nextflow.enable.dsl=2
workflow DECONV_INPUTS{
    // This is a function that prpeares the inputs for deconvolution that is utilised twice.
    take:
        cellbender_path
        prepare_inputs
    main:
        cellbender_path
            .map{row->tuple(row[0], file("${row[1]}".replaceFirst(/.*results/,"${params.outdir}")))}
            .set{ch_experiment_filth5} // this channel is used for task 'split_donor_h5ad'
            
        prepare_inputs.out.ch_experiment_bam_bai_barcodes.map { experiment, bam, bai, barcodes -> tuple(experiment,
                            bam,
                            bai)}.set{pre_ch_experiment_bam_bai_barcodes}
        cellbender_path.map{row->tuple(row[0], file("${row[1]}".replaceFirst(/.*results/,"${params.outdir}")+'/barcodes.tsv.gz'))}.set{barcodes}
        channel__file_paths_10x= cellbender_path.map{row->tuple(row[0],
                                                    file("${row[1]}".replaceFirst(/.*results/,"${params.outdir}")+'/barcodes.tsv.gz'),
                                                    file("${row[1]}".replaceFirst(/.*results/,"${params.outdir}")+'/features.tsv.gz'),
                                                    file("${row[1]}".replaceFirst(/.*results/,"${params.outdir}")+'/matrix.mtx.gz'))}
        pre_ch_experiment_bam_bai_barcodes.combine(barcodes, by: 0).set{ch_experiment_bam_bai_barcodes}
    emit:
        channel__file_paths_10x
        ch_experiment_bam_bai_barcodes
        ch_experiment_filth5
}