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
include { RETRIEVE_RECOURSES;RETRIEVE_RECOURSES_TEST_DATASET } from "$projectDir/modules/local/retrieve_recourses/retrieve_recourses"
include {RSYNC_RESULTS_REMOVE_WORK_DIR} from "$projectDir/modules/local/rsync_results_remove_work_dir/main"
include {celltype} from "$projectDir/subworkflows/celltype"
include {qc} from "$projectDir/subworkflows/qc"
include {dummy_filtered_channel} from "$projectDir/modules/local/merge_samples/functions"
include {CLUSTERING; CLUSTERING as CLUSTERING_HARMONY; CLUSTERING as CLUSTERING_BBKNN;} from "$projectDir/modules/local/clustering/main"
include { MATCH_GT_VIREO; GT_MATCH_POOL_IBD } from "$projectDir/modules/local/genotypes/main"
include { YASCP_INPUTS } from "$projectDir/modules/local/prepere_yascp_inputs/main"
include {MULTIPLET} from "$projectDir/subworkflows/doublet_detection"
include { match_genotypes } from "$projectDir/subworkflows/match_genotypes"
include { metadata_posthoc;replace_donors_posthoc } from "$projectDir/modules/local/report_update/main"
include {data_handover} from "$projectDir/subworkflows/data_handover"
include { CREATE_ARTIFICIAL_BAM_CHANNEL } from "$projectDir/modules/local/create_artificial_bam_channel/main"
include { TRANSFER;SUMMARY_STATISTICS_PLOTS } from "$projectDir/modules/local/summary_statistics_plots/main"
include {SUBSET_WORKF; JOIN_STUDIES_MERGE} from "$projectDir/modules/local/subset_genotype/main"
include {VIREO} from "$projectDir/modules/local/vireo/main"
include {capture_cellbender_files} from "$projectDir/modules/local/cellbender/functions"
include { DECONV_INPUTS } from "$projectDir/subworkflows/prepare_inputs"
include { prepare_inputs } from "$projectDir/subworkflows/prepare_inputs"
include { SPLIT_DONOR_H5AD } from "$projectDir/modules/local/split_donor_h5ad/main"
include {REPLACE_GT_DONOR_ID2 } from "$projectDir/modules/local/genotypes/main"
include {CAPTURE_VIREO } from "$projectDir/modules/local/vireo/main"
include {VIREO_GT_FIX_HEADER; VIREO_ADD_SAMPLE_PREFIX; MERGE_GENOTYPES_IN_ONE_VCF as MERGE_GENOTYPES_IN_ONE_VCF_INFERED; MERGE_GENOTYPES_IN_ONE_VCF as MERGE_GENOTYPES_IN_ONE_VCF_SUBSET} from "$projectDir/modules/local/genotypes/main"
include {ENHANCE_STATS_GT_MATCH } from "$projectDir/modules/local/genotypes/main"
include {collect_file} from "$projectDir/modules/local/collect_file/main"
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

workflow WORK_DIR_REMOVAL{
    // This process should be run with caution as it will remove the work directory and copy the results as an actual files


    out_ch = params.outdir
    ? Channel.fromPath(params.outdir, checkIfExists:true)
    : Channel.fromPath("${launchDir}/${outdir}")

    RSYNC_RESULTS_REMOVE_WORK_DIR(out_ch,params.tmpdir)
}

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
    qc(file__anndata_merged,file__cells_filtered,gt_outlier_input) //This runs the Clusterring and qc assessments of the datasets.    
}


workflow GT_MATCH{

    // 2) All the vcfs provided to us. 
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

workflow GENOTYPE_UPDATE{

    if (params.reference_assembly_fasta_dir=='https://yascp.cog.sanger.ac.uk/public/10x_reference_assembly'){
        RETRIEVE_RECOURSES()  
        genome = RETRIEVE_RECOURSES.out.reference_assembly
    }else{
        genome = "${params.reference_assembly_fasta_dir}"
    }

    // For Freeze1 we take the existing datasets and cp -as results folder so we can start from a breakpoint in pipeline
    // We rerun the GT match for all tranches as this has changed significantly since the beggining.
    myFileChannel = Channel.fromPath( "${params.outdir}/deconvolution/vireo_raw/*/GT_donors.vireo.vcf.gz" )
    myFileChannel.map{row -> tuple(row[-2], row)}.set{vireo_out_sample_donor_vcf}

    if (params.genotype_input.run_with_genotype_input) {
            // We have to produce a single vcf file for each individual pool.
            // Therefore we create 2 channels:
            // 1) All the expected vcf ids listed in the donor table
            Channel.fromPath(params.input_data_table,      
            followLinks: true,
            checkIfExists: true
            ).splitCsv(header: true, sep: '\t').map { row -> tuple(row.experiment_id, row.donor_vcf_ids) }
            .set { donors_in_pools }
            
            // 2) All the vcfs provided to us. 
            Channel.fromPath(
            params.genotype_input.tsv_donor_panel_vcfs,
            followLinks: true,
            checkIfExists: true
            ).splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) }
            .set { ch_ref_vcf }

            // This will subsequently result in a joint vcf file for all the cohorts listed for each of the pools that can be used in VIREO and/or GT matching algorythm.
            SUBSET_WORKF(ch_ref_vcf,donors_in_pools,'AllExpectedGT',genome)
            merged_expected_genotypes = SUBSET_WORKF.out.merged_expected_genotypes
            MERGE_GENOTYPES_IN_ONE_VCF_SUBSET(SUBSET_WORKF.out.study_merged_vcf.collect(),'subset')

        }else{
            ch_ref_vcf = Channel.of()
    }
    // tuple val(pool_id), path("${vireo_fixed_vcf}"), path("${vireo_fixed_vcf}.tbi"), emit: gt_pool
    // RERUN CellSNP
    input_channel = Channel.fromPath(params.input_data_table, followLinks: true, checkIfExists: true)
    CREATE_ARTIFICIAL_BAM_CHANNEL(input_channel)
    ch_experiment_bam_bai_barcodes = CREATE_ARTIFICIAL_BAM_CHANNEL.out.ch_experiment_bam_bai_barcodes
    ch_experiment_bam_bai_barcodes.view()
    channel_input_data_table = Channel.fromPath(params.input_data_table, followLinks: true, checkIfExists: true)
    channel_input_data_table
        .splitCsv(header: true, sep: params.input_tables_column_delimiter)
        .map{row->tuple(row.experiment_id, row.n_pooled)}
        .set{ch_experiment_npooled}    
    ch_experiment_bam_bai_barcodes.combine(ch_experiment_npooled, by: 0).set{cellsnp_with_npooled}
    // val(samplename), path(bam_file), path(bai_file), path(barcodes_tsv_gz),val(n_pooled)
    // [CRD_CMB13086620, /lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/fetch/Cardinal_45596_Aug_22_2022/nf_irods_to_lustre/pipelines/../../results/iget_study_cellranger/6776/CRD_CMB13086620/cellranger_CRD_CMB13086620/raw_feature_bc_matrix/../possorted_genome_bam.bam, /lustre/scratch123/hgi/mdt1/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/Cardinal_45596_Aug_22_2022/results_ct_adaptiveqc/results/deconvolution/vireo/CRD_CMB13086620/donor_ids.tsv, 12]
    if (params.existing_cellsnp != ''){
        log.info('Capturing some of the existing CELLSNP files')
        capture_cellsnp_files(params.existing_cellsnp)
        capture_cellsnp_files.out.cellsnp_loc.splitCsv(header: false, sep: ' ')
            .map{row->tuple(row[0], "${row[1]}")}
            .set{cellsnp_output_dir1}
        cellsnp_output_dir1.join(cellsnp_with_npooled, remainder: true).set{filter_channel}
        filter_channel.filter{ it[1] == null }.map{row -> tuple(row[0], row[2],row[3],row[4],row[5])}.set{cellsnp_with_npooled}
        cellsnp_with_npooled.view()
    }else{
        cellsnp_output_dir1 = Channel.of()
    }
    CELLSNP(cellsnp_with_npooled,
        Channel.fromPath(params.cellsnp.vcf_candidate_snps).collect())

    CAPTURE_VIREO(params.existing_vireo)
    CAPTURE_VIREO.out.vireo_loc.splitCsv(header: false, sep: ' ')
        .map{row->tuple(row[0], "${row[1]}", "${row[2]}")}
        .set{gt_pool}
    gt_pool
        .combine(ch_ref_vcf)
        .set { gt_math_pool_against_panel_input }
    
    vcf_input = Channel.fromPath(
                params.genotype_input.tsv_donor_panel_vcfs,
                followLinks: true,
                checkIfExists: true
            )

    vcf_input.splitCsv(header: true, sep: '\t')
            .map { row -> tuple(row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) }
            .set { ch_ref_vcf }
    match_genotypes(vireo_out_sample_donor_vcf,merged_expected_genotypes,gt_pool,gt_math_pool_against_panel_input,genome,ch_ref_vcf)
    ENHANCE_STATS_GT_MATCH(match_genotypes.out.donor_match_table_enhanced)
    collect_file(ENHANCE_STATS_GT_MATCH.out.assignments.collect(),"assignments_all_pools.tsv",params.outdir+'/deconvolution/vireo_processed',1,'')
    assignments_all_pools = collect_file.out.output_collection
    // We start the pipeline from the pre_qc breakpoint.
    emit:
        assignments_all_pools
}

workflow REPORT_UPDATE{
    // We use this entry point to update the reports upon running some individual processes that have already completed.
    input_channel = Channel.fromPath(params.input_data_table, followLinks: true, checkIfExists: true)
    CREATE_ARTIFICIAL_BAM_CHANNEL(input_channel)
    bam_split_channel = CREATE_ARTIFICIAL_BAM_CHANNEL.out.ch_experiment_bam_bai_barcodes
    ch_poolid_csv_donor_assignments = CREATE_ARTIFICIAL_BAM_CHANNEL.out.ch_poolid_csv_donor_assignments
    out_ch = params.outdir
    ? Channel.fromPath(params.outdir, checkIfExists:true)
    : Channel.fromPath("${launchDir}/${outdir}")
    // // 1) The pihat values were impemented posthoc, hence we are runing this on each of the independent tranches. 
    // myFileChannel = Channel.fromPath( "${params.outdir}/deconvolution/vireo/*/GT_donors.vireo.vcf.gz" )
    // myFileChannel.map{row -> tuple(row[-2], row)}.set{vireo_out_sample_donor_vcf}
    // match_genotypes(vireo_out_sample_donor_vcf)
    // match_genotypes.out.out_finish_val.set{o1}
    // // updating the Metadata if something new has been fetched,
    // // UKBB is sending us samples once a week and sometimes the sample mappings may be present at a later date, hence we update previously run samples accordingly.
    // Channel.from([["${params.RUN}",out_ch]]).set{update_input_channel}
    // // // We sometimes aslo chnge apporach in the data fetch and we need to add in some extra metadata
    // // metadata_posthoc(update_input_channel)
    // // metadata_posthoc.out.dummy_out.set{o3}
    // replace_donors_posthoc(update_input_channel)
    // replace_donors_posthoc.out.dummy_out.set{o2}
    // o1.mix(o2).last().set{o3}
    o3 = Channel.of('dummys')
    // Once everything is updated we need to make sure that the dataon the website and in the cardinal analysis foder is accurate and up to date, hence we rerun the data_handover scripts.
    // data_handover(params.outdir,
    //             process_finish_check_channel,
    //             ch_poolid_csv_donor_assignments,
    //             bam_split_channel) 



    data_handover(out_ch,
                input_channel,
                o3,
                ch_poolid_csv_donor_assignments,
                bam_split_channel) 
    // SUMMARY_STATISTICS_PLOTS(params.outdir,o3,params.input_data_table)
    // TRANSFER(SUMMARY_STATISTICS_PLOTS.out.summary_plots,params.rsync_to_web_file,params.outdir)
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

