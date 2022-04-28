
//include { prepare_inputs } from "$projectDir/subworkflows/prepare_inputs"
    //include { cellbender_deconv } from "$projectDir/subworkflows/cellbender"
include { main_deconvolution } from "$projectDir/subworkflows/main_deconvolution"
include { MERGE_SAMPLES } from "$projectDir/modules/nf-core/modules/merge_samples/main"

workflow deconvolution_module { 
    take:
        ch_experiment_bam_bai_barcodes
        prepare_inputs
        ch_experiment_filth5

    main:
        main_deconvolution(ch_experiment_bam_bai_barcodes, // activate this to run deconvolution pipeline
            prepare_inputs.out.ch_experiment_npooled,
            ch_experiment_filth5,
            prepare_inputs.out.ch_experiment_donorsvcf_donorslist,channel__file_paths_10x)

        MERGE_SAMPLES(main_deconvolution.out.out_h5ad,
                main_deconvolution.out.vireo_out_sample__exp_summary_tsv,
                'h5ad')

        file__anndata_merged = MERGE_SAMPLES.out.file__anndata_merged
        file__cells_filtered = MERGE_SAMPLES.out.file__cells_filtered

    emit:
        file__anndata_merged
        file__cells_filtered
}
