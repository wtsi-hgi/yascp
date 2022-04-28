
// Load base.config by default for all pipelines - typically included in the nextflow config.
include { CELLBENDER } from '../modules/nf-core/modules/cellbender/main'
include { DECONV_INPUTS } from "${projectDir}/subworkflows/prepare_inputs/deconvolution_inputs"

workflow cellbender {
    take:
        ch_experimentid_paths10x_raw
		ch_experimentid_paths10x_filtered
        channel__metadata
    main:
        log.info params.input_data_table
        log.info """---Running Cellbender pipeline ---"""
        CELLBENDER(ch_experimentid_paths10x_raw,ch_experimentid_paths10x_filtered,channel__metadata)
        results_list = CELLBENDER.out.results_list
        cellbender_path=CELLBENDER.out.cellbender_path
    emit:
        results_list
        cellbender_path
}

workflow cellbender_deconvolution {
    take:
        input
        prepare_inputs
    main:
           // // Removing the background using cellbender which is then used in the deconvolution.
        if (input == 'cellbender'){
            log.info ' ---- using cellbender to remove background---'
            cellbender(prepare_inputs.out.ch_experimentid_paths10x_raw,
                prepare_inputs.out.ch_experimentid_paths10x_filtered,
                prepare_inputs.out.channel__metadata)
            log.info ' ---- Out results - cellbender to remove background---'
            DECONV_INPUTS(cellbender.out.cellbender_path,prepare_inputs)

            channel__file_paths_10x = DECONV_INPUTS.out.channel__file_paths_10x
            ch_experiment_bam_bai_barcodes= DECONV_INPUTS.out.ch_experiment_bam_bai_barcodes
            ch_experiment_filth5= DECONV_INPUTS.out.ch_experiment_filth5

        }else if (input == 'existing_cellbender'){
            // Here we are using the existing cellbender from a different run, Nothe that the structure of the cellbender folder should be same as produced by this pipeline.
            log.info ' ---- using existing cellbender output for deconvolution---'
            capture_cellbender_files(params.cellbender_location,"${params.output_dir}/nf-preprocessing")
            DECONV_INPUTS(capture_cellbender_files.out.celbender_path,prepare_inputs)

            channel__file_paths_10x = DECONV_INPUTS.out.channel__file_paths_10x
            ch_experiment_bam_bai_barcodes= DECONV_INPUTS.out.ch_experiment_bam_bai_barcodes
            ch_experiment_filth5= DECONV_INPUTS.out.ch_experiment_filth5
        }
        else if (input == 'cellranger'){
            // This is where we skip the cellbender and use the cellranger filtered datasets.
            log.info '--- using cellranger filtered data instead of cellbender (skipping cellbender)---'
            ch_experiment_bam_bai_barcodes=prepare_inputs.out.ch_experiment_bam_bai_barcodes
            ch_experiment_filth5=prepare_inputs.out.ch_experiment_filth5
            channel__file_paths_10x=prepare_inputs.out.channel__file_paths_10x
        }
        else{
            log.info '--- input mode is not selected - please choose --- (existing_cellbender | cellbender | cellranger)'
        }
    emit:
        channel__file_paths_10x 
        ch_experiment_bam_bai_barcodes
        ch_experiment_filth5
}