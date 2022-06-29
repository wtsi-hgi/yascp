#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//include { MATCH_GT_VIREO } from '../modules/nf-core/modules/genotypes/main.nf'
include { match_genotypes } from '../subworkflows/match_genotypes.nf'
params.vireo.run_gtmatch_aposteriori = true


workflow TEST_MATCH_GT_VIREO {
  ch_pool_id = Channel.value( 'Card_Val12156644' )
  ch_vireo_vcf = Channel.fromPath(
    '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/Cherry_ELGH/results/deconvolution/vireo/Card_Val12156644/GT_donors.vireo.vcf.gz'
    )
  ch_pool_id.combine(ch_vireo_vcf)
  .set{ ch_pool_id_vireo_vcf }
  ch_pool_id_vireo_vcf.subscribe { println "TEST_MATCH_GT_VIREO: ch_pool_id_vireo_vcf = ${it}" }

  Channel.fromPath(
    //'/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/Cherry_ELGH/gtmatch/1000G_full_GRCh38.srt.vcf.gz'
    '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/genotypes/filtered_genotypes/vcf/merged_bcf/GT_AF_ELGH_Concat.bcf.gz'
  )
  .map { fnam -> tuple(fnam, "${fnam}.csi") }
  .set { ch_ref_vcf }

  ch_ref_vcf.subscribe { println "TEST_MATCH_GT_VIREO: ch_ref_vcf = ${it}" }
  match_genotypes(ch_pool_id_vireo_vcf, ch_ref_vcf)
}
