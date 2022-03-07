#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { MATCH_GT_VIREO } from '../modules/nf-core/modules/genotypes/main.nf'

workflow TEST_MATCH_GT_VIREO {
  ch_pool_id = Channel.value( 'Card_Val12156644' )
  ch_vireo_vcf = Channel.fromPath( '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/Cherry_ELGH/results/deconvolution/vireo/Card_Val12156644/GT_donors.vireo.vcf.gz' )
  ch_pool_id.merge(ch_vireo_vcf).set(ch_pool_id_vireo_vcf)
  Channel.fromPath(
  '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/Cherry_ELGH/gtmatch/1000G_full_GRCh38.srt.vcf.gz'
  )
  .map { file -> tuple(file, "${file}.tbi")}
  .subscribe { println "TEST_MATCH_GT_VIREO: ${it}" }
  .set { ch_ref_vcf }

  MATCH_GT_VIREO(ch_pool_id_vireo_vcf, ch_ref_vcf)
}
