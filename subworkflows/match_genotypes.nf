// match deconvoluted donors by genotype to a reference panel

include { MATCH_GT_VIREO; GT_MATCH_POOL_IBD } from '../modules/nf-core/modules/genotypes/main'

workflow match_genotypes {
  take:
    ch_pool_id_vireo_vcf

  main:
    Channel.fromPath(
      params.genotype_input.tsv_donor_panel_vcfs,
      followLinks: true,
      checkIfExists: true
    )
    .splitCsv(header: true, sep: '\t')
    .map { row -> tuple(row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) }
    .set { ch_ref_vcf }

    MATCH_GT_VIREO(ch_pool_id_vireo_vcf, ch_ref_vcf)

    // compare genotypes within a pool (identity by descent)
    GT_MATCH_POOL_IBD(ch_pool_id_vireo_vcf)

  emit:
    pool_id_donor_assignments_csv = MATCH_GT_VIREO.out.pool_id_donor_assignments_csv
    donor_match_table = MATCH_GT_VIREO.out.donor_match_table
}
