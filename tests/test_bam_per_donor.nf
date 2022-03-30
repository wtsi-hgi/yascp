#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

//include { SPLIT_CELL_BARCODES_PER_DONOR } from '../modules/local/cellranger_bam_per_donor.nf'
include { split_bam_by_donor } from '../modules/local/cellranger_bam_per_donor.nf'

workflow TEST_SPLIT_BAM_PER_DONOR {
  Channel.fromPath(
    '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/hp3_dev/elgh_yascp/testdata/testdata.tsv'
    ).splitCsv(header: true, sep: '\t')
    .map{ row->tuple("${row.pool_id}", file("${row.bam_file}"), file("${row.donor_ids}")) }
    .set{ ch_sample_possorted_bam_vireo_donor_ids }

    ch_sample_possorted_bam_vireo_donor_ids
      .subscribe onNext: { println "ch_sample_possorted_bam_vireo_donor_ids = ${it}" }, onComplete: { println "Done.\n"}

  //SPLIT_CELL_BARCODES_PER_DONOR(ch_sample_possorted_bam_vireo_donor_ids)

  split_bam_by_donor(ch_sample_possorted_bam_vireo_donor_ids)
}
