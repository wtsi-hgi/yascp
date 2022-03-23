
// Load base.config by default for all pipelines - typically included in the nextflow config.
include { main_deconvolution } from './main_deconvolution'
include { MATCH_GT_VIREO } from '../modules/nf-core/modules/genotypes/main'
include { SPLIT_BAM_PER_DONOR } from '../modules/local/cellranger_bam_per_donor'

workflow deconvolution {
    take:
        ch_experiment_bam_bai_barcodes
        ch_experiment_npooled
        ch_experiment_filth5
        ch_experiment_donorsvcf_donorslist
        channel__file_paths_10x
    main:
        // // run main deconvolution pipeline on prepared input channels:
        main_deconvolution(ch_experiment_bam_bai_barcodes,
                ch_experiment_npooled,
                ch_experiment_filth5,
                ch_experiment_donorsvcf_donorslist,channel__file_paths_10x)
        out_h5ad = main_deconvolution.out.out_h5ad
        vireo_out_sample__exp_summary_tsv=main_deconvolution.out.vireo_out_sample__exp_summary_tsv
        vireo_out_sample_donor_vcf = main_deconvolution.out.vireo_out_sample_donor_vcf
    emit:
        out_h5ad
        vireo_out_sample__exp_summary_tsv
        vireo_out_sample_donor_vcf
        sample_possorted_bam_vireo_donor_ids = main_deconvolution.out.sample_possorted_bam_vireo_donor_ids
}

workflow match_genotypes {
  take:
    ch_pool_id_vireo_vcf
    ch_ref_vcf

  main:
    MATCH_GT_VIREO(ch_pool_id_vireo_vcf, ch_ref_vcf)

  emit:
    donor_assignments_csv = MATCH_GT_VIREO.out.donor_match_table
    gtcheck_results_txt = MATCH_GT_VIREO.out.gtcheck_out
}

workflow split_bam_by_donor {
  take:
    sample_possorted_bam_vireo_donor_ids

  main:
    SPLIT_BAM_PER_DONOR(sample_possorted_bam_vireo_donor_ids)
}
