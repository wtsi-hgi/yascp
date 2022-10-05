#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/yascp
========================================================================================
    Github : https://github.com/nf-core/yascp
    Website: https://nf-co.re/yascp
    Slack  : https://nfcore.slack.com/channels/yascp
----------------------------------------------------------------------------------------
*/


nextflow.enable.dsl = 2

/*
========================================================================================
    VALIDATE & PRINT PARAMETER SUMMARY
========================================================================================
*/


/*
========================================================================================
  NAMED WORKFLOWS FOR TESTS
========================================================================================
*/

//include { TEST_MATCH_GENOTYPES } from './tests/test_genotypes'
//include { TEST_GT_ASSIGN } from './tests/test_gt_assign'
//include { TEST_GTCHECK } from './tests/test_gtcheck'
//include { TEST_MATCH_GT_VIREO } from './tests/test_genotypes'
//include { TEST_SPLIT_BAM_PER_DONOR } from './tests/test_bam_per_donor'
//include { TEST_ENCRYPT_DIR} from './tests/test_encryption'
include { TEST_SUBSET_GENOTYPES } from './tests/test_subset_genotypes'
include { match_genotypes } from './subworkflows/match_genotypes'

workflow NF_CORE_TEST {
  //println "**** running NF_CORE_TEST::TEST_MATCH_GENOTYPES"
  //TEST_MATCH_GENOTYPES()
  //println "**** running NF_CORE_TEST::TEST_GT_ASSIGN"
  //TEST_GT_ASSIGN()
  //println "**** running NF_CORE_TEST::TEST_GTCHECK"
  //TEST_GTCHECK()
  //println "**** running NF_CORE_TEST::TEST_MATCH_GT_VIREO"
  //TEST_MATCH_GT_VIREO()
  //println "**** running NF_CORE_TEST::TEST_SPLIT_BAM_PER_DONOR"
  //TEST_SPLIT_BAM_PER_DONOR()
  //println "**** running NF_CORE_TEST::TEST_ENCRYPT_DIR"
  //TEST_ENCRYPT_DIR()
  println "**** running NF_CORE_TEST::TEST_SUBSET_GENOTYPES"
  TEST_SUBSET_GENOTYPES()
}


workflow REPORT_UPDATE{
    // [CRD_CMB13102396, /lustre/scratch123/hgi/mdt1/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/Cardinal_45717_Sep_03_2022/work/28/9c0b54dbae1243523421cf09c2a40f/vireo_CRD_CMB13102396/GT_donors.vireo.vcf.gz]
    // match_genotypes(vireo_out_sample_donor_vcf)
    // qc/Cardinal_45717_Sep_03_2022/results/deconvolution/vireo/CRD_CMB13102395/GT_donors.vireo.vcf.gz
    log.info "#### LETS run the Genotype matching algorythm and PiHat Calculations #####"
    myFileChannel = Channel.fromPath( "${params.outdir}/deconvolution/vireo/*/GT_donors.vireo.vcf.gz" )
    myFileChannel.map{row -> tuple(row[-2], row)}.set{vireo_out_sample_donor_vcf}
    match_genotypes(vireo_out_sample_donor_vcf)
    
}


/*
========================================================================================
    NAMED WORKFLOW FOR PIPELINE
========================================================================================
*/

include { SCDECON } from './workflows/yascp'

//
// WORKFLOW: Run main nf-core/yascp analysis pipeline
//
workflow NFCORE_SCDECON {
    SCDECON ()
}

/*
========================================================================================
    RUN ALL WORKFLOWS
========================================================================================
*/

//
// WORKFLOW: Execute a single named workflow for the pipeline
// See: https://github.com/nf-core/rnaseq/issues/619
//
workflow {
    NFCORE_SCDECON ()
}

/*
========================================================================================
    THE END
========================================================================================
*/
