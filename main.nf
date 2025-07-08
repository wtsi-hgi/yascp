#!/usr/bin/env nextflow
/*
========================================================================================
    wtsi-hgi/yascp
========================================================================================
    Github : https://github.com/wtsi-hgi/yascp
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

include { YASCP } from "$projectDir/workflows/yascp"
include { RETRIEVE_RECOURSES; RETRIEVE_RECOURSES_TEST_DATASET } from "$projectDir/modules/local/retrieve_recourses/retrieve_recourses"
include { RSYNC_RESULTS_REMOVE_WORK_DIR} from "$projectDir/modules/local/rsync_results_remove_work_dir/main"
include { celltype} from "$projectDir/subworkflows/celltype"
include { qc_and_integration } from "$projectDir/subworkflows/qc_and_integration"
include { dummy_filtered_channel } from "$projectDir/modules/local/merge_samples/functions"
include { CLUSTERING; CLUSTERING as CLUSTERING_HARMONY; 
          CLUSTERING as CLUSTERING_BBKNN} from "$projectDir/modules/local/clustering/main"
include { MATCH_GT_VIREO; GT_MATCH_POOL_IBD } from "$projectDir/modules/local/genotypes/main"
include { YASCP_INPUTS } from "$projectDir/modules/local/prepere_yascp_inputs/main"
include { MULTIPLET } from "$projectDir/subworkflows/doublet_detection"
include { match_genotypes } from "$projectDir/subworkflows/match_genotypes"
include { metadata_posthoc;
          replace_donors_posthoc } from "$projectDir/modules/local/report_update/main"
include { data_handover } from "$projectDir/subworkflows/data_handover"
include { CREATE_ARTIFICIAL_BAM_CHANNEL } from "$projectDir/modules/local/create_artificial_bam_channel/main"
include { TRANSFER;
          SUMMARY_STATISTICS_PLOTS } from "$projectDir/modules/local/summary_statistics_plots/main"
include { SUBSET_WORKF; 
          JOIN_STUDIES_MERGE} from "$projectDir/modules/local/subset_genotype/main"
include { VIREO } from "$projectDir/modules/local/vireo/main"
include {capture_cellbender_files} from "$projectDir/modules/local/cellbender/functions"
include { DECONV_INPUTS } from "$projectDir/subworkflows/prepare_inputs"
include { prepare_inputs } from "$projectDir/subworkflows/prepare_inputs"
include { SPLIT_DONOR_H5AD } from "$projectDir/modules/local/split_donor_h5ad/main"
include { REPLACE_GT_DONOR_ID2 } from "$projectDir/modules/local/genotypes/main"
include { CAPTURE_VIREO } from "$projectDir/modules/local/vireo/main"
include { VIREO_GT_FIX_HEADER; VIREO_ADD_SAMPLE_PREFIX; 
        MERGE_GENOTYPES_IN_ONE_VCF as MERGE_GENOTYPES_IN_ONE_VCF_INFERED; 
        MERGE_GENOTYPES_IN_ONE_VCF as MERGE_GENOTYPES_IN_ONE_VCF_SUBSET } from "$projectDir/modules/local/genotypes/main"
include { ENHANCE_STATS_GT_MATCH } from "$projectDir/modules/local/genotypes/main"
include { collect_file} from "$projectDir/modules/local/collect_file/main"
include { CELLSNP;capture_cellsnp_files } from "$projectDir/modules/local/cellsnp/main"
include { CONVERT_H5AD_TO_MTX } from "$projectDir/modules/local/convert_h5ad_to_mtx/main"

////// WORKFLOW: Run main wtsi-hgi/yascp analysis pipeline
// This is the default entry point, we have others to update ceirtain parts of the results. 
// Please go to ./workflows/yascp to see the main Yascp workflow.

workflow MAIN {

    out_ch = params.outdir
            ? Channel.fromPath(params.outdir, checkIfExists:true)
            : Channel.fromPath("${launchDir}/${outdir}")
    if (params.profile=='test_full'){
        RETRIEVE_RECOURSES_TEST_DATASET(out_ch)
        input_channel = RETRIEVE_RECOURSES_TEST_DATASET.out.input_channel
        vcf_inputs = RETRIEVE_RECOURSES_TEST_DATASET.out.vcf_inputs
        vcf_inputs.splitCsv(header: true, sep: '\t')
                    .map { row -> tuple(row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) }
                    .set { vcf_inputs }
    }else{
        input_channel = Channel.fromPath(params.input_data_table, followLinks: true, checkIfExists: true)
        if (params.genotype_input.run_with_genotype_input) {
            vcf_inputs = Channel.fromPath(
                params.genotype_input.tsv_donor_panel_vcfs,
                followLinks: true,
                checkIfExists: true
            )
        vcf_inputs.splitCsv(header: true, sep: '\t')
                    .map { row -> tuple(row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) }
                    .set { vcf_inputs }
        }else{
            vcf_inputs = Channel.of()
        }
    }

    input_channel.collectFile(name: "$params.outdir/yascp_inputs.tsv")
    YASCP ('default',input_channel,vcf_inputs)

}

workflow {
    MAIN ()
}

////////////////////////////
// END OF MAIN ENTRANCE IN WORKFLOWS - No need to look beyond ulness you are looking at a specific entry point.
//////////////////////////////

workflow JUST_CELLTYPES{
    file__anndata_merged = Channel.from(params.file__anndata_merged)
    celltype(file__anndata_merged,'celltype_mode')
}

workflow JUST_CELLBENDER{
    // here we are skipping everything downstram and are only performing cellbender opperations
    params.do_deconvolution = false
    params.celltype_assignment.run_celltype_assignment = false
    params.skip_qc = true
    params.skip_handover = true
    params.skip_merge = true
    MAIN ()
}


workflow JUST_DOUBLETS{
    file__anndata_merged = Channel.from(params.file__anndata_merged)
    MULTIPLET(
        file__anndata_merged,'doublet_mode'
    )
}

workflow JUST_RECLUSTER{
    // This workflow entry point is for integrations only. Users can provide a h5ad file that can be clustered. 
    file__anndata_merged = Channel.from(params.file__anndata_merged)
    gt_outlier_input = Channel.from("$projectDir/assets/fake_file.fq")
    file__anndata_merged.subscribe { println "file__anndata_merged: $it" }
    dummy_filtered_channel(file__anndata_merged,params.id_in)
    file__cells_filtered = dummy_filtered_channel.out.anndata_metadata
    qc_and_integration(file__anndata_merged,file__cells_filtered,gt_outlier_input) //This runs the Clusterring and qc assessments of the datasets.    
}


workflow GT_MATCH{
    Channel.fromPath(
    params.genotype_input.tsv_donor_panel_vcfs,
    followLinks: true,
    checkIfExists: true
    ).splitCsv(header: true, sep: '\t')
    .map { row -> tuple(params.file_name,file("${params.vcf}"),file("${params.vcf}.tbi"),row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) }
    .set { gt_math_pool_against_panel_input }

    MATCH_GT_VIREO(gt_math_pool_against_panel_input)
    ENHANCE_STATS_GT_MATCH(MATCH_GT_VIREO.out.donor_match_table_with_pool_id,params.input_data_table)
}


workflow WORK_DIR_REMOVAL{
    // This process should be run with caution as it will remove the work directory and copy the results as an actual files
    out_ch = params.outdir
    ? Channel.fromPath(params.outdir, checkIfExists:true)
    : Channel.fromPath("${launchDir}/${outdir}")
    RSYNC_RESULTS_REMOVE_WORK_DIR(out_ch,params.tmpdir)
}


workflow.onComplete{
    log.info "Pipeline completed at: $workflow.complete"
    log.info "Command line: $workflow.commandLine"
    log.info "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
    if (workflow.success){
        println "Execution successful"
        if (params.remove_work_dir){
            println "lets remove ${params.tmpdir}"
            log.info "You have selected \"remove_work_dir = true\"; will therefore remove work dirs of all tasks"
            def proc = "bash ${projectDir}/bin/del_work_dirs.sh ${params.tmpdir}".execute()
            def b = new StringBuffer()
            proc.consumeProcessErrorStream(b)
            log.info proc.text
            log.info b.toString() 
        }
    }
}

