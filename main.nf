#!/usr/bin/env nextflow
/*
========================================================================================
    wtsi-hgi/yascp
========================================================================================
    Github : https://github.com/wtsi-hgi/yascp
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

// --- MODULE/SUBWORKFLOW IMPORTS ---
include { YASCP } from "$projectDir/workflows/yascp"
include { RETRIEVE_RECOURSES; RETRIEVE_RECOURSES_TEST_DATASET } from "$projectDir/modules/local/retrieve_resources/retrieve_resources"
include { RSYNC_RESULTS_REMOVE_WORK_DIR} from "$projectDir/modules/local/rsync_results_remove_work_dir/main"
include { CELLTYPE} from "$projectDir/subworkflows/celltype"
include { QC_AND_INTEGRATION } from "$projectDir/subworkflows/qc_and_integration"
include { DUMMY_FILTERED_CHANNEL } from "$projectDir/modules/local/merge_samples/functions"
include { MATCH_GT_VIREO } from "$projectDir/modules/local/genotypes/main"
include { MULTIPLET } from "$projectDir/subworkflows/doublet_detection"
include { ENHANCE_STATS_GT_MATCH } from "$projectDir/modules/local/genotypes/main"

/**
 * Helper function to parse VCF input TSVs consistently across workflows
 */
def parseVcfTsv(tsv_path) {
    return Channel.fromPath(tsv_path, followLinks: true, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row -> 
            tuple(row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) 
        }
}

/*
========================================================================================
    MAIN WORKFLOW ENTRY
========================================================================================
*/

workflow MAIN {
    // Determine Output Directory Channel
    def outdir_path = params.outdir.value ?: "${launchDir}/${outdir}"
    out_ch = Channel.fromPath(outdir_path, checkIfExists: true)

    // Input Handling: Test vs. Production
    if (params.profile == 'test_full') {
        RETRIEVE_RECOURSES_TEST_DATASET(out_ch)
        input_channel = RETRIEVE_RECOURSES_TEST_DATASET.out.input_channel
        vcf_inputs = RETRIEVE_RECOURSES_TEST_DATASET.out.vcf_inputs
            .splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) }
    } else {
        input_channel = Channel.fromPath(params.input_data_table.value, followLinks: true, checkIfExists: true)
        
        vcf_inputs = params.genotype_input.run_with_genotype_input
            ? parseVcfTsv(params.genotype_input.tsv_donor_panel_vcfs)
            : Channel.empty()
    }

    // Log tracking
    input_channel.collectFile(name: "${params.outdir.value}/yascp_inputs.tsv")
    
    // Execute Core Workflow
    YASCP('default', input_channel, vcf_inputs)
}

// Default implicit workflow
workflow {
    MAIN()
}

/*
========================================================================================
    SPECIALIZED ENTRY POINTS (Sub-Workflows)
========================================================================================
*/

workflow JUST_CELLTYPES {
    CELLTYPE(Channel.fromPath(params.file__anndata_merged.value), 'celltype_mode')
}

workflow JUST_CELLBENDER {
    // Override params for lean execution
    params.do_deconvolution.value = false
    params.celltype_assignment.run_celltype_assignment.value = false
    params.skip_qc.value = true
    params.skip_handover.value = true
    //params.skip_merge = true
    MAIN()
}

workflow JUST_DOUBLETS {
    MULTIPLET(Channel.fromPath(params.file__anndata_merged.value), 'doublet_mode')
}

workflow JUST_RECLUSTER {
    ch_merged = Channel.fromPath(params.file__anndata_merged.value)
    
    // Optimization: Using empty list instead of "fake_file" assets where applicable
    gt_outlier_input = Channel.value([]) 
    
    DUMMY_FILTERED_CHANNEL(ch_merged, params.id_in.value)
    QC_AND_INTEGRATION(ch_merged, DUMMY_FILTERED_CHANNEL.out.anndata_metadata, gt_outlier_input)
}

workflow GT_MATCH {
    // Standardizing GT matching logic
    parseVcfTsv(params.genotype_input.tsv_donor_panel_vcfs)
        .map { label, vcf, csi -> 
            tuple(params.file_name.value, file(params.vcf.value), file("${params.vcf.value}.tbi"), label, vcf, csi)
        }
        .set { gt_match_input }

    MATCH_GT_VIREO(gt_match_input)
    ENHANCE_STATS_GT_MATCH(MATCH_GT_VIREO.out.donor_match_table_with_pool_id, params.input_data_table.value)
}

workflow WORK_DIR_REMOVAL {
    def outdir_path = params.outdir.value ?: "${launchDir}/${outdir}"
    RSYNC_RESULTS_REMOVE_WORK_DIR(Channel.fromPath(outdir_path, checkIfExists: true), params.tmpdir.value)
}

/*
========================================================================================
    COMPLETION HANDLER
========================================================================================
*/

workflow.onComplete {
    log.info """
    Pipeline execution summary
    ---------------------------
    Completed at: ${workflow.complete}
    Duration    : ${workflow.duration}
    Success     : ${workflow.success}
    workDir     : ${workflow.workDir}
    exit status : ${workflow.exitStatus}
    """
    
    if (workflow.success && params.remove_work_dir.value) {
        log.info "Cleaning up work directory: ${params.tmpdir.value}"
        "bash ${projectDir}/bin/del_work_dirs.sh ${params.tmpdir.value}".execute().waitFor()
    }
}