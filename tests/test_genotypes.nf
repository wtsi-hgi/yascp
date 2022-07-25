#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { match_genotypes } from '../subworkflows/match_genotypes.nf'
params.vireo.run_gtmatch_aposteriori = true

workflow TEST_MATCH_GENOTYPES
{
  ch_pool_id_vireo_vcf = Channel.from(
  ['CRD_CMB12979963', '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results_copy/deconvolution/vireo/CRD_CMB12979963/GT_donors.vireo.vcf.gz'],
  ['CRD_CMB12979969', '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results_copy/deconvolution/vireo/CRD_CMB12979969/GT_donors.vireo.vcf.gz'],
  ['CRD_CMB12979966', '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results_copy/deconvolution/vireo/CRD_CMB12979966/GT_donors.vireo.vcf.gz'],
  ['CRD_CMB12979972', '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results_copy/deconvolution/vireo/CRD_CMB12979972/GT_donors.vireo.vcf.gz'],
  ['CRD_CMB12979964', '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results_copy/deconvolution/vireo/CRD_CMB12979964/GT_donors.vireo.vcf.gz'],
  ['CRD_CMB12979965', '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results_copy/deconvolution/vireo/CRD_CMB12979965/GT_donors.vireo.vcf.gz'],
  ['CRD_CMB12979968', '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results_copy/deconvolution/vireo/CRD_CMB12979968/GT_donors.vireo.vcf.gz'],
  ['CRD_CMB12979971', '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results_copy/deconvolution/vireo/CRD_CMB12979971/GT_donors.vireo.vcf.gz'],
  ['CRD_CMB12979967', '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/UKBB_ELGH_5th_July_2022/results_copy/deconvolution/vireo/CRD_CMB12979967/GT_donors.vireo.vcf.gz']
  ).map { nam, filpath -> [nam, file(filpath)] }

  ch_pool_id_vireo_vcf.subscribe { println "TEST_MATCG_GENOTYPE: ch_pool_id_vireo_vcf = ${it}" }

  match_genotypes(ch_pool_id_vireo_vcf)
  match_genotypes.out.pool_id_donor_assignments_csv.view()
}
