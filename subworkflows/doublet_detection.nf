include { SCRUBLET } from "$projectDir/modules/local/scrublet/main"
include { DOUBLET_DETECTION } from "$projectDir/modules/local/doubletdetection/main"
include { DOUBLET_DECON } from "$projectDir/modules/local/doubletdecon/main"
include { SC_DBLFINDER } from "$projectDir/modules/local/scDblFinder/main"
include { DOUBLET_FINDER} from "$projectDir/modules/local/doubletfinder/main"
include { SCDS } from "$projectDir/modules/local/scds/main"
include { SPLIT_CITESEQ_GEX; 
          SPLIT_CITESEQ_GEX as SPLIT_CITESEQ_GEX_FILTERED } from "$projectDir/modules/local/citeseq/main"
include { CONVERT_MTX_TO_H5AD; CONVERT_H5AD_TO_MTX } from "$projectDir/modules/local/convert_h5ad_to_mtx/main"
include { SPLIT_BATCH_H5AD } from "$projectDir/modules/local/split_batch_h5ad/main"


process MERGE_DOUBLET_RESULTS{
    tag "${experiment_id}"
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "${params.yascp_container}"

    } else {
        container "wtsihgi/nf_scrna_qc:6bb6af5"
    }

    publishDir path: "${params.outdir}/doublet_detection/droplet_type_distribution", 
               mode: "${params.copy_mode}", 
               pattern: "*.png",
               overwrite: "true"
    publishDir path: "${params.outdir}/doublet_detection/doublet_results_combined", 
               mode: "${params.copy_mode}", 
               pattern: "*.tsv",
               overwrite: "true"

    input:
        tuple(
            val(experiment_id),
            path(all_doublet_files)
        )

    output:
        path("*.png") optional true
        tuple val(experiment_id), path("${experiment_id}__doublet_results_combined.tsv"), emit: result

    script:
        
        """
            echo 'Lets combine'
            combine_doublets.py --list ${all_doublet_files}
            ln -s all_doublet_results_combined.tsv ${experiment_id}__doublet_results_combined.tsv
            doublet_plots.py
            mv droplet_type_distribution.png ${experiment_id}__droplet_type_distribution.png
        """
}

workflow MULTIPLET {
    take:
        channel__file_paths_10x
        mode
    main:
        // Identify multiplets.
        log.info '---Identifying doublets and multiplets---'
        input_channel = Channel.of()
        if (mode=='yascp_full'){
            gex_h5ad = CONVERT_MTX_TO_H5AD(channel__file_paths_10x).gex_h5ad
        }else{
            log.info '---Splitting the assignment for each batch---'
            SPLIT_BATCH_H5AD(channel__file_paths_10x,params.doublet_celltype_split_column)
            SPLIT_BATCH_H5AD.out.sample_file
                .splitCsv(header: true, sep: "\t", by: 1)
                .map{row -> tuple(row.experiment_id, file(row.h5ad_filepath))}.set{gex_h5ad}      
            SPLIT_BATCH_H5AD.out.sample_file
                .splitCsv(header: true, sep: "\t", by: 1)
                .map{row -> file(row.h5ad_filepath)}.set{gex_h5ad_file} 
            channel__file_paths_10x_pre = CONVERT_H5AD_TO_MTX(gex_h5ad_file).channel__file_paths_10x    
            channel__file_paths_10x =  channel__file_paths_10x_pre
                .map{row -> tuple(
                row[0],
                file("${row[1]}/barcodes.tsv.gz"),
                file("${row[1]}/features.tsv.gz"),
                file("${row[1]}/matrix.mtx.gz")
            )} 
        }

        if (params.filter_multiplets.scrublet.run_process){
            SCRUBLET(
                channel__file_paths_10x,
                params.filter_multiplets.expected_multiplet_rate,
                params.filter_multiplets.scrublet.n_simulated_multiplet,
                params.filter_multiplets.scrublet.multiplet_threshold_method,
                params.filter_multiplets.scrublet.scale_log10
            )
            input_channel = input_channel.mix(SCRUBLET.out.result)
        }else{
            out = Channel.of()
            input_channel = Channel.of()
        }

        if (params.filter_multiplets.doubletDetection.run_process){
            DOUBLET_DETECTION(channel__file_paths_10x) //Done
            input_channel = input_channel.mix(DOUBLET_DETECTION.out.result)
        }else{
            out = Channel.of()
        }

        if (params.filter_multiplets.doubletDecon.run_process){
            DOUBLET_DECON(gex_h5ad) //Done
            input_channel = input_channel.mix(DOUBLET_DECON.out.result)
        }else{
            out = Channel.of()
        }

        if (params.filter_multiplets.scDblFinder.run_process){
            SC_DBLFINDER(channel__file_paths_10x)
            input_channel = input_channel.mix(SC_DBLFINDER.out.result)
        }else{
            out = Channel.of()
        }

        if (params.filter_multiplets.scds.run_process){
            SCDS(channel__file_paths_10x)
            input_channel = input_channel.mix(SCDS.out.result)
        }else{
            out = Channel.of()
        }

        if (params.filter_multiplets.doubletFinder.run_process){
            DOUBLET_FINDER(gex_h5ad,params.filter_multiplets.expected_multiplet_rate) //Done
            input_channel = input_channel.mix(DOUBLET_FINDER.out.result)
            
        }else{
            out = Channel.of()
        }

        input_channel2 = input_channel.groupTuple()
        MERGE_DOUBLET_RESULTS(input_channel2) 


    emit:
        scrublet_paths = MERGE_DOUBLET_RESULTS.out.result
}

