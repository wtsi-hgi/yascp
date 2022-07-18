// match deconvoluted donors by genotype to a reference panel

include { MATCH_GT_VIREO } from '../modules/nf-core/modules/genotypes/main'

workflow match_genotypes {
  take:
    ch_pool_id_vireo_vcf

  main:
    Channel.fromPath(
      params.tsv_donor_panel_vcfs,
      followLinks: true,
      checkIfExists: true
    )
    .splitCsv(header: true, sep: '\t')
    .map { row -> tuple(row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) }
    .set { ch_ref_vcf }

    match_gt_vireo(ch_pool_id_vireo_vcf, ch_ref_vcf)

  emit:
    pool_id_donor_assignments_csv = MATCH_GT_VIREO.out.pool_id_donor_assignments_csv
}
