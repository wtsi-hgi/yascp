/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/

ch_multiqc_config        = file("$projectDir/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config) : Channel.empty()

/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

def modules = params.modules.clone()

include { GET_SOFTWARE_VERSIONS } from "$projectDir/modules/local/get_software_versions" addParams( options: [publish_files : ['tsv':'']] )
include { main_deconvolution } from "$projectDir/subworkflows/main_deconvolution"
include { match_genotypes } from "$projectDir/subworkflows/match_genotypes"
include {cellbender} from "$projectDir/subworkflows/cellbender"
include {qc} from "$projectDir/subworkflows/qc"
include {data_handover} from "$projectDir/subworkflows/data_handover"
include { prepare_inputs } from "$projectDir/subworkflows/prepare_inputs"
include { DECONV_INPUTS } from "$projectDir/subworkflows/prepare_inputs/deconvolution_inputs"
include {MERGE_SAMPLES} from "$projectDir/modules/nf-core/modules/merge_samples/main"
include {MULTIPLET} from "../modules/nf-core/modules/multiplet/main"
include {dummy_filtered_channel} from "../modules/nf-core/modules/merge_samples/functions"
include {capture_cellbender_files} from "../modules/nf-core/modules/cellbender/functions"

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SCDECON {

    // ###################################
    // ################################### Readme
    // Step1. CELLBENDER
    // There are 3 modes of running YASCP pipeline: 
    // (option 1) users can run it from 10x data and use cellbender -  params.input == 'cellbender' 
    // (option 2) users can run it from existing cellbender if the analysis has already been performed -  params.input == 'existing_cellbender' : note a specific folder structure is required
    // (option 3) users can run it from cellranger - skipping the cellbender. params.input == 'cellranger'
    // ###################################
    // ###################################

    if (!params.skip_preprocessing.value){
        input_channel = Channel.fromPath(params.input_data_table, followLinks: true, checkIfExists: true)
        // prepearing the inputs from a standard 10x dataset folders.
        prepare_inputs(input_channel)
        log.info 'The preprocessing has been already performed, skipping directly to h5ad input'
           // // Removing the background using cellbender which is then used in the deconvolution.
        if (params.input == 'cellbender'){
            log.info ' ---- using cellbender to remove background---'
            cellbender(prepare_inputs.out.ch_experimentid_paths10x_raw,
                prepare_inputs.out.ch_experimentid_paths10x_filtered,prepare_inputs.out.channel__metadata)
            log.info ' ---- Out results - cellbender to remove background---'
            DECONV_INPUTS(cellbender.out.cellbender_path,prepare_inputs)
            channel__file_paths_10x = DECONV_INPUTS.out.channel__file_paths_10x
            ch_experiment_bam_bai_barcodes= DECONV_INPUTS.out.ch_experiment_bam_bai_barcodes
            ch_experiment_filth5= DECONV_INPUTS.out.ch_experiment_filth5
        }else if (params.input == 'existing_cellbender'){
            // Here we are using the existing cellbender from a different run, Nothe that the structure of the cellbender folder should be same as produced by this pipeline.
            log.info ' ---- using existing cellbender output for deconvolution---'
            capture_cellbender_files(params.cellbender_location,"${params.output_dir}/nf-preprocessing")
            DECONV_INPUTS(capture_cellbender_files.out.celbender_path,prepare_inputs)
            channel__file_paths_10x = DECONV_INPUTS.out.channel__file_paths_10x
            ch_experiment_bam_bai_barcodes= DECONV_INPUTS.out.ch_experiment_bam_bai_barcodes
            ch_experiment_filth5= DECONV_INPUTS.out.ch_experiment_filth5
        }
        else if (params.input == 'cellranger'){
            // This is where we skip the cellbender and use the cellranger filtered datasets.
            log.info '--- using cellranger filtered data instead of cellbender (skipping cellbender)---'
            ch_experiment_bam_bai_barcodes=prepare_inputs.out.ch_experiment_bam_bai_barcodes
            ch_experiment_filth5=prepare_inputs.out.ch_experiment_filth5
            channel__file_paths_10x=prepare_inputs.out.channel__file_paths_10x
        }
        else{
            log.info '--- input mode is not selected - please choose --- (existing_cellbender| cellbender | cellranger)'
        }
        // ###################################
        // ################################### Readme
        // Step2. DECONVOLUTION
        // When thepreprocessing with cellbender or cellranger is finalised then we can do the deconvolution of samples. This can also be skipped if the samples are not multiplexed.
        // However if the number of individuals is specified as 1 the deconvolution withh be skipped anyways, but we will apply scrubblet to remove dublicates.
        // Suggestion is to still run deconvolution so that dublicates are removed.
        // ###################################
        // ###################################

        if (params.do_deconvolution){
            main_deconvolution(ch_experiment_bam_bai_barcodes, // activate this to run deconvolution pipeline
                prepare_inputs.out.ch_experiment_npooled,
                ch_experiment_filth5,
                prepare_inputs.out.ch_experiment_donorsvcf_donorslist,channel__file_paths_10x)
            MERGE_SAMPLES(main_deconvolution.out.out_h5ad,main_deconvolution.out.vireo_out_sample__exp_summary_tsv,'h5ad')
            
        }else{
            channel__metadata = prepare_inputs.out.channel__metadata
            MERGE_SAMPLES(channel__file_paths_10x,channel__metadata,'barcodes')
        }
        // TODO: Here add a fundtion to take an extra h5ad and merge it together with the current run. This will be required for the downstream analysis when we want to integrate multiple datasets
        file__anndata_merged = MERGE_SAMPLES.out.file__anndata_merged
        file__cells_filtered = MERGE_SAMPLES.out.file__cells_filtered
    }else{
        // This option skips all the deconvolution and and takes a preprocessed yascp h5ad file to run the downstream clustering and celltype annotation.
        log.info '''----Skipping Preprocessing since we already have prepeared h5ad input file----'''
        file__anndata_merged = Channel.from(params.skip_preprocessing.file__anndata_merged)
        if (params.skip_preprocessing.file__cells_filtered ==''){
            log.info '''--- No cells filtered input ----'''
            dummy_filtered_channel(file__anndata_merged,params.skip_preprocessing.id_in)
            file__cells_filtered = dummy_filtered_channel.out.anndata_metadata
        }else{
            file__cells_filtered = Channel.from(params.skip_preprocessing.file__cells_filtered)
        }
    }

    // ###################################
    // ################################### Readme
    // Step3. QC METRICS, CELLTYPE ASSIGNMENT and CLUSTERIN
    // After background removal and demultiplexing we perform qc metrics and clustering of the processed cells. 
    // This step of the pipeline also performs celltype assignments and removes cells that fail adaptive filtering.
    // ###################################
    // ###################################

    // qc(file__anndata_merged,file__cells_filtered) //This runs the Clusterring and qc assessments of the datasets.

    // The idea is to also run eQTL analysis, however this is currently not implemented as part of this pipeline.
    // // // Performing eQTL mapping.

    // ###################################
    // ################################### Readme
    // Step4. SUMMARY STATISTICS, DONOR SPLITTING and WEB TRANSFER
    // Once all the processes are done we gather the data in a summary folder in results. This is mainly done for the purposes of reporting on the web and o split the donor based metrics. 
    // ###################################
    // ###################################

    // data_handover("${workDir}/../${params.output_dir}",qc.out.LI) //This part gathers the plots for the reporting in a Summary folder. If run through gitlab CI it will triger the data transfer to web.
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
