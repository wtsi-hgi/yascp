/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)



// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist

// Check mandatory parameters
// if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input samplesheet not specified!' }

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

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )
include {deconvolution} from "$projectDir/subworkflows/deconvolution"
include {cellbender} from "$projectDir/subworkflows/cellbender"
// include {qc} from "$projectDir/subworkflows/qc"
include { prepare_inputs } from "$projectDir/subworkflows/prepare_inputs"
//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { INPUT_CHECK } from '../subworkflows/local/input_check' addParams( options: [:] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

def multiqc_options   = modules['multiqc']
multiqc_options.args += params.multiqc_title ? Utils.joinModuleArgs(["--title \"$params.multiqc_title\""]) : ''

//
// MODULE: Installed directly from nf-core/modules
//
// include { FASTQC  } from '../modules/nf-core/modules/fastqc/main'  addParams( options: modules['fastqc'] )
// include { MULTIQC } from '../modules/nf-core/modules/multiqc/main' addParams( options: multiqc_options   )

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/

// Info required for completion email and summary
def multiqc_report = []

workflow SCDECON {

    // prepearing the inputs from a standard 10x dataset folders.
	prepare_inputs(Channel.fromPath(params.input_data_table, followLinks: true, checkIfExists: true))

    // // Removing the background using cellbender which is then used in the deconvolution.
    if (params.input == 'cellbender'){
        log.info ' ---- using cellbender to remove background---'
        cellbender(prepare_inputs.out.ch_experimentid_paths10x_raw,
            prepare_inputs.out.ch_experimentid_paths10x_filtered)
        log.info ' ---- Out results - cellbender to remove background---'
        cellbender.out.results_list.view() //TODO
        cellbender.out.results_list
            .map{experiment, path -> tuple(experiment, path+'/cellbender_FPR_0pt01_filtered.h5')}
            .set{ch_experiment_filth5} // this channel is used for task 'split_donor_h5ad'

        prepare_inputs.out.ch_experiment_bam_bai_barcodes.map { experiment, bam, bai, barcodes -> tuple(experiment,
                    bam,
                    bai)}.set{pre_ch_experiment_bam_bai_barcodes}

        cellbender.out.results_list
            .map{experiment, path -> tuple(experiment, file(path+'/cellbender-FPR_0pt01-filtered_10x_mtx/barcodes.tsv.gz'))}.set{barcodes}

        pre_ch_experiment_bam_bai_barcodes.combine(barcodes, by: 0).set{ch_experiment_bam_bai_barcodes}
        ch_experiment_bam_bai_barcodes.view()
        // TODO - have to alter 2 channels - ch_experiment_filth5 and ch_experiment_bam_bai_barcodes based on the cellbender output paths cellbender.out.results_list


    }else if (params.input == 'existing_cellbender'){
        log.info ' ---- using existing cellbender output for deconvolution---'
        Channel.fromPath(params.cellbender_file, followLinks: true, checkIfExists: true)
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.experiment_id, row.data_path_10x_format)}
            .map{experiment, path -> tuple(experiment, path.replaceFirst(/_10x_mtx/,".h5"))}
            .map{experiment, path -> tuple(experiment, path.replaceFirst(/cellbender-FPR/,"cellbender_FPR"))}
            .map{experiment, path -> tuple(experiment, path.replaceFirst(/-filtered.h5/,"_filtered.h5"))}.set{ch_experiment_filth5} // this channel is used for task 'split_donor_h5ad'

        prepare_inputs.out.ch_experiment_bam_bai_barcodes.map { experiment, bam, bai, barcodes -> tuple(experiment,
							bam,
							bai)}.set{pre_ch_experiment_bam_bai_barcodes}


        Channel.fromPath(params.cellbender_file, followLinks: true, checkIfExists: true)
            .splitCsv(header: true, sep: params.input_tables_column_delimiter).map{row->tuple(row.experiment_id, file(row.data_path_10x_format+'/barcodes.tsv.gz'))}.set{barcodes}

        pre_ch_experiment_bam_bai_barcodes.combine(barcodes, by: 0).set{ch_experiment_bam_bai_barcodes}

    }else{
        log.info '--- using cellranger filtered data instead of cellbender (skipping cellbender)---'
        ch_experiment_bam_bai_barcodes=prepare_inputs.out.ch_experiment_bam_bai_barcodes
        ch_experiment_filth5=prepare_inputs.out.ch_experiment_filth5
    }

    // Performing Deconvolution of Samples.
    deconvolution(ch_experiment_bam_bai_barcodes, // activate this to run deconvolution pipeline
		       prepare_inputs.out.ch_experiment_npooled,
		       ch_experiment_filth5,
		       prepare_inputs.out.ch_experiment_donorsvcf_donorslist)


    // Performing QC metrics -
    // TODO we may want to split this in the Clustring, QC, Celltype, Web_transfer
    // qc()



    // Performing eQTL mapping.
    // This part will contain code from Hannes and the potentially additional LIMIX runs.






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
