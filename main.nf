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
include { TEST_SUBSET_GENOTYPES } from "$projectDir/tests/test_subset_genotypes"
include { match_genotypes } from "$projectDir/subworkflows/match_genotypes"
include { metadata_posthoc;replace_donors_posthoc } from "$projectDir/modules/local/report_update/main"
include {data_handover} from "$projectDir/subworkflows/data_handover"
include { CREATE_ARTIFICIAL_BAM_CHANNEL } from "$projectDir/modules/local/create_artificial_bam_channel/main"
include { TRANSFER;SUMMARY_STATISTICS_PLOTS } from "$projectDir/modules/nf-core/modules/summary_statistics_plots/main"

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
    // We use this entry point to update the reports upon running some individual processes that have already completed.
    input_channel = Channel.fromPath(params.input_data_table, followLinks: true, checkIfExists: true)
    CREATE_ARTIFICIAL_BAM_CHANNEL(input_channel)
    bam_split_channel = CREATE_ARTIFICIAL_BAM_CHANNEL.out.ch_experiment_bam_bai_barcodes
    ch_poolid_csv_donor_assignments = CREATE_ARTIFICIAL_BAM_CHANNEL.out.ch_poolid_csv_donor_assignments
    

    // 1) The pihat values were impemented posthoc, hence we are runing this on each of the independent tranches. 
    myFileChannel = Channel.fromPath( "${params.outdir}/deconvolution/vireo/*/GT_donors.vireo.vcf.gz" )
    myFileChannel.map{row -> tuple(row[-2], row)}.set{vireo_out_sample_donor_vcf}
    match_genotypes(vireo_out_sample_donor_vcf)

    // updating the Metadata if something new has been fetched,
    // UKBB is sending us samples once a week and sometimes the sample mappings may be present at a later date, hence we update previously run samples accordingly.
    Channel.from([["${params.RUN}","${params.output_dir}"]]).set{update_input_channel}
    // We sometimes aslo chnge apporach in the data fetch and we need to add in some extra metadata
    // metadata_posthoc(update_input_channel)
    replace_donors_posthoc(update_input_channel)
    process_finish_check_channel = replace_donors_posthoc.out.dummy_out
    // Once everything is updated we need to make sure that the dataon the website and in the cardinal analysis foder is accurate and up to date, hence we rerun the data_handover scripts.
    // data_handover(params.output_dir,
    //             process_finish_check_channel,
    //             ch_poolid_csv_donor_assignments,
    //             bam_split_channel) 

    SUMMARY_STATISTICS_PLOTS(params.output_dir,process_finish_check_channel,params.input_data_table)
    TRANSFER(SUMMARY_STATISTICS_PLOTS.out.summary_plots)
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
