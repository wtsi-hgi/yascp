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
