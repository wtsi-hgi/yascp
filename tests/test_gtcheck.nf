#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
  VIREO_GT_FIX_HEADER;
  GT_MATCH_POOL_AGAINST_PANEL;
  ASSIGN_DONOR_FROM_PANEL;
  ASSIGN_DONOR_OVERALL
  } from '../modules/nf-core/modules/genotypes/main.nf'
include { match_genotypes } from '../subworkflows/match_genotypes.nf'
//params.vireo.run_gtmatch_aposteriori = true


workflow TEST_GTCHECK {
//ch_pool_id = Channel.value( 'CRD_CMB12979963', 'CRD_CMB12979969' )
//ch_vireo_vcf = Channel.fromPath(
//  '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results/deconvolution/vireo/CRD_CMB12979963/GT_donors.vireo.vcf.gz',
//  '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results/deconvolution/vireo/CRD_CMB12979969/GT_donors.vireo.vcf.gz'
  //'/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/hp3_dev/elgh_220526_mutref/results/deconvolution/vireo/CRD_CMB12922402/GT_donors.vireo.vcf.gz'
//)
//  ch_pool_id.combine(ch_vireo_vcf)
//  .set{ ch_pool_id_vireo_vcf }

  ch_pool_id_vireo_vcf = Channel.from(
    ['CRD_CMB12979963', '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results/deconvolution/vireo/CRD_CMB12979963/GT_donors.vireo.vcf.gz'],
    ['CRD_CMB12979969', '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results/deconvolution/vireo/CRD_CMB12979969/GT_donors.vireo.vcf.gz'],
    ['CRD_CMB12979966', '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results/deconvolution/vireo/CRD_CMB12979966/GT_donors.vireo.vcf.gz'],
    ['CRD_CMB12979972', '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results/deconvolution/vireo/CRD_CMB12979972/GT_donors.vireo.vcf.gz'],
    ['CRD_CMB12979964', '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results/deconvolution/vireo/CRD_CMB12979964/GT_donors.vireo.vcf.gz'],
    ['CRD_CMB12979965', '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results/deconvolution/vireo/CRD_CMB12979965/GT_donors.vireo.vcf.gz'],
    ['CRD_CMB12979968', '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results/deconvolution/vireo/CRD_CMB12979968/GT_donors.vireo.vcf.gz'],
    ['CRD_CMB12979971', '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results/deconvolution/vireo/CRD_CMB12979971/GT_donors.vireo.vcf.gz'],
    ['CRD_CMB12979967', '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results/deconvolution/vireo/CRD_CMB12979967/GT_donors.vireo.vcf.gz']
    ).map { nam, filpath -> [nam, file(filpath)] }

  ch_pool_id_vireo_vcf.subscribe { println "TEST_GT_CHECK: ch_pool_id_vireo_vcf = ${it}" }


  Channel.fromPath(
    '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/hp3_dev/elgh_yascp/donor_panels_vcf_paths.tsv',
    followLinks: true,
    checkIfExists: true
  )
  .splitCsv(header: true, sep: '\t')
  .map { row -> tuple(row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) }
  .set { ch_ref_vcf }

  ch_ref_vcf.subscribe { println "TEST_GTCHECK: ch_ref_vcf = ${it}" }


  VIREO_GT_FIX_HEADER(ch_pool_id_vireo_vcf)

  ch_gt_pool = VIREO_GT_FIX_HEADER.out.gt_pool
  ch_gt_pool
    .combine(ch_ref_vcf)
    .set { ch_gt_pool_ref_vcf }

  ch_gt_pool_ref_vcf.subscribe { println "GT_MATCH_POOL_AGAINST_PANEL: ch_gt_pool_ref_vcf = ${it}\n" }
  //ch_gt_pool_ref_vcf.view()
  GT_MATCH_POOL_AGAINST_PANEL(ch_gt_pool_ref_vcf)
  // group by panel id
  GT_MATCH_POOL_AGAINST_PANEL.out.gtcheck_results
    .groupTuple()
    .set { gt_check_by_panel }

  gt_check_by_panel.subscribe { println "ASSIGN_DONOR_FROM_PANEL: gt_check_by_panel = ${it}\n" }
  ASSIGN_DONOR_FROM_PANEL(gt_check_by_panel)

  ASSIGN_DONOR_FROM_PANEL.out.gtcheck_assignments
    .groupTuple()
    .set{ ch_donor_assign_panel }

  ch_donor_assign_panel.subscribe {println "ASSIGN_DONOR_OVERALL: ch_donor_assign_panel = ${it}\n"}
  ASSIGN_DONOR_OVERALL(ch_donor_assign_panel)

  ASSIGN_DONOR_OVERALL.out.donor_assignments.view()
  //match_genotypes(ch_pool_id_vireo_vcf, ch_ref_vcf)
}
