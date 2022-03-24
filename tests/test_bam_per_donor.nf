#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//include { SPLIT_CELL_BARCODES_PER_DONOR } from '../modules/local/cellranger_bam_per_donor.nf'
include { split_bam_by_donor } from '../modules/local/cellranger_bam_per_donor.nf'

workflow TEST_SPLIT_BAM_PER_DONOR {
  Channel.fromPath(
    '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/hp3_dev/elgh_yascp/testdata/ELGH_VAL11509202/cellranger_ELGH_VAL11509202/possorted_genome_bam.bam'
    ).set{ ch_cellranger_possorted_bam }
  Channel.fromPath(
    '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/hp3_dev/elgh_yascp/testdata/ELGH_VAL11509202/donor_ids.tsv'
    ).set { ch_vireo_donor_barcode_tsv }

  Channel.value( 'ELGH_VAL11509202' )
    .combine(ch_cellranger_possorted_bam)
    .combine(ch_vireo_donor_barcode_tsv)
    .set {ch_sample_possorted_bam_vireo_donor_ids}

    ch_sample_possorted_bam_vireo_donor_ids
      .subscribe { println "ch_sample_possorted_bam_vireo_donor_ids = ${it}" }

  //SPLIT_CELL_BARCODES_PER_DONOR(ch_sample_possorted_bam_vireo_donor_ids)
  split_bam_by_donor(ch_sample_possorted_bam_vireo_donor_ids)
}
