#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { SPLIT_BAM_PER_DONOR } from '../modules/local/cellranger_bam_per_donor.nf'


workflow TEST_SPLIT_BAM_PER_DONOR {
  Channel.fromPath(
    '/lustre/scratch123/hgi/mdt1/projects/ukbb_scrna/pipelines/Pilot_UKB/fetch/ELGH/results/iget_study_cellranger/6634/ELGH_VAL11509202/cellranger_ELGH_VAL11509202/possorted_genome_bam.bam'
    ).set{ ch_cellranger_possorted_bam }
  Channel.fromPath(
    '/lustre/scratch123/hgi/mdt1/projects/ukbb_scrna/pipelines/hp3_dev/elgh_yascp/results/deconvolution/vireo/ELGH_VAL11509202/donor_ids.tsv'
    ).set { ch_vireo_donor_barcode_tsv }
  SPLIT_BAM_PER_DONOR(ch_vireo_donor_barcode_tsv, ch_cellranger_possorted_bam)
}
