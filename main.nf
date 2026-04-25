#!/usr/bin/env nextflow
/*
========================================================================================
    nf-core/yascp
========================================================================================
    Github : https://github.com/nf-core/yascp
    Slack  : https://nfcore.slack.com/channels/yascp
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2
include { YASCP } from "$projectDir/workflows/yascp"
include { RETRIEVE_RECOURSES;RETRIEVE_RECOURSES_TEST_DATASET } from "$projectDir/subworkflows/local/retrieve_recourses"
include {RSYNC_RESULTS_REMOVE_WORK_DIR} from "$projectDir/modules/local/rsync_results_remove_work_dir/main"
include {celltype} from "$projectDir/subworkflows/celltype"
include {qc} from "$projectDir/subworkflows/qc"
include {dummy_filtered_channel} from "$projectDir/modules/nf-core/modules/merge_samples/functions"
include {CLUSTERING; CLUSTERING as CLUSTERING_HARMONY; CLUSTERING as CLUSTERING_BBKNN;} from "$projectDir/modules/nf-core/modules/clustering/main"
include { MATCH_GT_VIREO; GT_MATCH_POOL_IBD } from "$projectDir/modules/nf-core/modules/genotypes/main"
include { YASCP_INPUTS } from "$projectDir/modules/nf-core/modules/prepere_yascp_inputs/main"
include {MULTIPLET} from "$projectDir/subworkflows/doublet_detection"
include { match_genotypes } from "$projectDir/subworkflows/match_genotypes"
include { metadata_posthoc;replace_donors_posthoc } from "$projectDir/modules/local/report_update/main"
include {data_handover} from "$projectDir/subworkflows/data_handover"
include { CREATE_ARTIFICIAL_BAM_CHANNEL } from "$projectDir/modules/local/create_artificial_bam_channel/main"
include { TRANSFER;SUMMARY_STATISTICS_PLOTS } from "$projectDir/modules/nf-core/modules/summary_statistics_plots/main"
include {SUBSET_WORKF; JOIN_STUDIES_MERGE} from "$projectDir/modules/nf-core/modules/subset_genotype/main"
include {VIREO} from "$projectDir/modules/nf-core/modules/vireo/main"
include {capture_cellbender_files} from "$projectDir/modules/nf-core/modules/cellbender/functions"
include { DECONV_INPUTS } from "$projectDir/subworkflows/prepare_inputs/deconvolution_inputs"
include { prepare_inputs } from "$projectDir/subworkflows/prepare_inputs"
include { SPLIT_DONOR_H5AD } from "$projectDir/modules/nf-core/modules/split_donor_h5ad/main"
include {REPLACE_GT_DONOR_ID2 } from "$projectDir/modules/nf-core/modules/genotypes/main"
include {CAPTURE_VIREO } from "$projectDir/modules/nf-core/modules/vireo/main"
include {VIREO_GT_FIX_HEADER; VIREO_ADD_SAMPLE_PREFIX; MERGE_GENOTYPES_IN_ONE_VCF as MERGE_GENOTYPES_IN_ONE_VCF_INFERED; MERGE_GENOTYPES_IN_ONE_VCF as MERGE_GENOTYPES_IN_ONE_VCF_SUBSET} from "$projectDir/modules/nf-core/modules/genotypes/main"
include {ENHANCE_STATS_GT_MATCH } from "$projectDir/modules/nf-core/modules/genotypes/main"
include {collect_file} from "$projectDir/modules/nf-core/modules/collect_file/main"
include { CELLSNP;capture_cellsnp_files } from "$projectDir/modules/nf-core/modules/cellsnp/main"
include { CONVERT_H5AD_TO_MTX } from "$projectDir/modules/local/convert_h5ad_to_mtx/main"


// Channel.of(params.outdir).mkdirs()

////// WORKFLOW: Run main nf-core/yascp analysis pipeline
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
    celltype(file__anndata_merged)
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
    // input_channel = Channel.fromPath(params.input_data_table, followLinks: true, checkIfExists: true)
    // YASCP_INPUTS(input_channel)
    // channel_input_data_table = YASCP_INPUTS.out.input_file_corectly_formatted


    // channel__file_paths_10x =  channel_input_data_table
    //     .splitCsv(header: true, sep: params.input_tables_column_delimiter)
    //     .map{row -> tuple(
    //     row.experiment_id,
    //     file("${row.data_path_10x_format}/filtered_feature_bc_matrix")
    // )}
    file__anndata_merged = Channel.from(params.file__anndata_merged)
    channel__file_paths_10x = CONVERT_H5AD_TO_MTX(file__anndata_merged).channel__file_paths_10x

    channel__file_paths_10x =  channel__file_paths_10x
        .map{row -> tuple(
        row[0],
        file("${row[1]}/barcodes.tsv.gz"),
        file("${row[1]}/features.tsv.gz"),
        file("${row[1]}/matrix.mtx.gz")
    )}

    MULTIPLET(
        channel__file_paths_10x
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


////// You do not need to concern about the workflows bellow as these are Cardinal Specific and used for development
/*
========================================================================================
  Below we have other workflow that are a versions of the Yascp to avoid ceirtain modules and update/validate the datasets
========================================================================================
*/



workflow FREEZE1_GENERATION{
    GENOTYPE_UPDATE()

    input_channel = Channel.fromPath(params.input_data_table, followLinks: true, checkIfExists: true)
    if (params.genotype_input.run_with_genotype_input) {
        vcf_inputs = Channel.fromPath(
            params.genotype_input.tsv_donor_panel_vcfs,
            followLinks: true,
            checkIfExists: true
        )
    }else{
        vcf_inputs = Channel.of()
    }
    vcf_input.splitCsv(header: true, sep: '\t')
                .map { row -> tuple(row.label, file(row.vcf_file_path), file("${row.vcf_file_path}.csi")) }
                .set { vcf_input }
    YASCP (GENOTYPE_UPDATE.out.assignments_all_pools,input_channel,vcf_inputs)



}


workflow TEST_CATCHE_ISSUES{
    GENOTYPE_UPDATE()
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


workflow TEST {
  //println "**** running TEST::TEST_MATCH_GENOTYPES"
  //TEST_MATCH_GENOTYPES()
  //println "**** running TEST::TEST_GT_ASSIGN"
  //TEST_GT_ASSIGN()
  //println "**** running TEST::TEST_GTCHECK"
  //TEST_GTCHECK()
  //println "**** running TEST::TEST_MATCH_GT_VIREO"
  //TEST_MATCH_GT_VIREO()
  //println "**** running TEST::TEST_SPLIT_BAM_PER_DONOR"
  //TEST_SPLIT_BAM_PER_DONOR()
  //println "**** running TEST::TEST_ENCRYPT_DIR"
  //TEST_ENCRYPT_DIR()
          if (params.reference_assembly_fasta_dir=='https://yascp.cog.sanger.ac.uk/public/10x_reference_assembly'){
            RETRIEVE_RECOURSES()  
            genome = RETRIEVE_RECOURSES.out.reference_assembly
        }else{
            genome = "${params.reference_assembly_fasta_dir}"
        }
        
    println "**** running TEST::TEST_SUBSET_GENOTYPES"
    //   TEST_SUBSET_GENOTYPES()
    input_channel = Channel.fromPath(params.input_data_table, followLinks: true, checkIfExists: true)
    
    prepare_inputs(input_channel)

    if(params.cellbender_location==''){
        cellbender_location = "${params.outdir}/dummy"
    }

    capture_cellbender_files(cellbender_location,"${params.outdir}/nf-preprocessing")
    DECONV_INPUTS(capture_cellbender_files.out.celbender_path,prepare_inputs)
    channel__file_paths_10x = DECONV_INPUTS.out.channel__file_paths_10x
    ch_experiment_filth5= DECONV_INPUTS.out.ch_experiment_filth5
    MULTIPLET(
        params.outdir,
        channel__file_paths_10x,
        params.sample_qc.cell_filters.filter_multiplets.expected_multiplet_rate,
        params.sample_qc.cell_filters.filter_multiplets.n_simulated_multiplet,
        params.sample_qc.cell_filters.filter_multiplets.multiplet_threshold_method,
        params.sample_qc.cell_filters.filter_multiplets.scale_log10
    )

    log.info "#### running DECONVOLUTION workflow #####"
    if (params.run_with_genotype_input) {
        // if (params.genotype_input.subset_genotypes){
            log.info "---We are subsetting genotypes----"
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
            SUBSET_WORKF(ch_ref_vcf,donors_in_pools)
            merged_expected_genotypes = SUBSET_WORKF.out.merged_expected_genotypes
            // merged_expected_genotypes.view()
            // SUBSET_GENOTYPE(ch_experiment_donorsvcf_donorslist.map { experiment, donorsvcf,donortbi, donorslist -> tuple(experiment,
            //                 file(donorsvcf),file(donortbi),
            //                 donorslist)})
        // }
    }

    channel_input_data_table = Channel.fromPath(params.input_data_table, followLinks: true, checkIfExists: true)
        
    if (params.existing_cellsnp != ''){
        expl1 = Channel.fromPath( "${params.existing_cellsnp}/*/*.vcf.gz")
        expl1.map{row -> tuple("${row[-2]}".replaceAll('cellsnp_',''), "${row}".replaceAll('/cellSNP.base.vcf.gz',''))}.set{cellsnp_output_dir}
        
    }else{
        log.inf('Running CELLSNP')
        channel_input_data_table
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.experiment_id, "${row.data_path_10x_format}/possorted_genome_bam.bam" ,row.data_path_10x_format+'/filtered_feature_bc_matrix/barcodes.tsv.gz')}
            .set{pre_ch_experiment_bam_barcodes}
        pre_ch_experiment_bam_barcodes
            .map { a,b,c -> tuple(a, file(b), file("${b}.bai"), file(c))}
            .set {ch_experiment_bam_bai_barcodes}

        CELLSNP(ch_experiment_bam_bai_barcodes,
            Channel.fromPath(params.cellsnp.vcf_candidate_snps).collect())
        cellsnp_output_dir = CELLSNP.out.cellsnp_output_dir
    }

    channel_input_data_table
        .splitCsv(header: true, sep: params.input_tables_column_delimiter)
        .map{row->tuple(row.experiment_id, row.n_pooled)}
        .set{ch_experiment_npooled}

    if (params.run_with_genotype_input) {
        cellsnp_output_dir.combine(ch_experiment_npooled, by: 0)
            .combine(merged_expected_genotypes, by: 0).set{full_vcf}  
    }else{
        log.info "-----running Vireo without genotype input----"
        cellsnp_output_dir.combine(ch_experiment_npooled, by: 0).set{full_vcf}
        full_vcf.map {experiment, cellsnp, npooled -> tuple(experiment, cellsnp, npooled,[],[])}.set{full_vcf}
    }

    full_vcf.filter { experiment, cellsnp, npooled, t,ti -> npooled != '1' }.set{full_vcf2}
    full_vcf.filter { experiment, cellsnp, npooled, t,ti -> npooled == '1' }.set{not_deconvoluted}
            
    VIREO(full_vcf2)  
    // To make things constent we still use donor0 ... donor1 input names, hence the folowing module to change this back to a default approach
    
    REPLACE_GT_DONOR_ID2(VIREO.out.all_required_data)

    vireo_out_sample_donor_vcf = REPLACE_GT_DONOR_ID2.out.infered_vcf
    vireo_out_sample_summary_tsv = REPLACE_GT_DONOR_ID2.out.sample_summary_tsv
    vireo_out_sample__exp_summary_tsv = REPLACE_GT_DONOR_ID2.out.sample__exp_summary_tsv
    vireo_out_sample_donor_ids = REPLACE_GT_DONOR_ID2.out.sample_donor_ids

    // vireo_out_sample_donor_vcf2 = VIREO.out.infered_vcf
    // vireo_out_sample_summary_tsv2 = VIREO.out.sample_summary_tsv
    // vireo_out_sample__exp_summary_tsv2 = VIREO.out.sample__exp_summary_tsv
    // vireo_out_sample_donor_ids2 = VIREO.out.sample_donor_ids

    not_deconvoluted.map{ experiment, donorsvcf, npooled,t,t2 -> tuple(experiment, 'None')}.set{not_deconvoluted2}
    file_cellmetadata = MULTIPLET.out.file__cellmetadata
    scrublet_paths = MULTIPLET.out.scrublet_paths
    split_channel = vireo_out_sample_donor_ids.combine(ch_experiment_filth5, by: 0)
    split_channel2 = not_deconvoluted2.combine(ch_experiment_filth5, by: 0)
    // combining these 2 channels in one
    split_channel3 = split_channel.mix(split_channel2)
    // adding the scrublet paths to the channel.
    split_channel4 = split_channel3.combine(scrublet_paths, by: 0)

    split_channel5 = split_channel4.map{
        val_sample, val_donor_ids_tsv, val_filtered_matrix_h5, path_scrublet ->
        [  val_sample,file(val_donor_ids_tsv),file(val_filtered_matrix_h5),path_scrublet,params.outdir]
    }
    
    SPLIT_DONOR_H5AD(split_channel5)
    if (params.run_with_genotype_input) {
        match_genotypes(vireo_out_sample_donor_vcf,merged_expected_genotypes)
    }

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

