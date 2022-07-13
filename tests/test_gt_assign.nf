#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {
  ASSIGN_DONOR_FROM_PANEL;
  ASSIGN_DONOR_OVERALL
  } from '../modules/nf-core/modules/genotypes/main.nf'

//params.vireo.run_gtmatch_aposteriori = true
params.fofn = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/hp3_dev/elgh_yascp/test/test_CRD_CMB12979963/fofn.tsv"

workflow TEST_GT_ASSIGN {

  ch_fofn_in = Channel.fromPath(params.fofn, followLinks: true, checkIfExists: true)
  ch_fofn_in.splitCsv(header: true, sep: '\t')
  .map{row->tuple(row.pool_panel_id, file(row.file_path))}
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
