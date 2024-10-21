include { SCRUBLET } from "$projectDir/modules/nf-core/modules/scrublet/main"
include { DOUBLET_DETECTION } from "$projectDir/modules/nf-core/modules/doubletdetection/main"
include { DOUBLET_DECON} from "$projectDir/modules/nf-core/modules/doubletdecon/main"
include {SC_DBLFINDER} from "$projectDir/modules/nf-core/modules/scDblFinder/main"
include { DOUBLET_FINDER} from "$projectDir/modules/nf-core/modules/doubletfinder/main"
include { SCDS} from "$projectDir/modules/nf-core/modules/scds/main"
include { SPLIT_CITESEQ_GEX; SPLIT_CITESEQ_GEX as SPLIT_CITESEQ_GEX_FILTERED } from "$projectDir/modules/nf-core/modules/citeseq/main"
include { CONVERT_MTX_TO_H5AD } from "$projectDir/modules/local/convert_h5ad_to_mtx/main"

def random_hex(n) {
    Long.toUnsignedString(new Random().nextLong(), n).toUpperCase()
}

process make_cellmetadata_pipeline_input {
    // Makes a input tsv file for the main pipeline.
    // ------------------------------------------------------------------------
    //cache false        // cache results from run

    publishDir  path: "${params.outdir}/multiplet.method=scrublet",
                saveAs: {filename -> filename.replaceAll("-", "")},
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        path("*multiplet_calls_published.txt")

    output:
        path('file_cellmetadata.tsv', emit: file__cellmetadata)

    script:

        """
        # Note: the default paste delim is tab
        cat *multiplet_calls_published.txt \
            | awk 'BEGIN{print "experiment_id\tdata_path_cellmetadata"}1' \
            > file_cellmetadata.tsv
        """
}


process MERGE_DOUBLET_RESULTS{
    tag "${experiment_id}"
    label 'process_medium'
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://yascp.cog.sanger.ac.uk/public/singularity_images/nf_scrna_qc_v3.img"
    } else {
        container "mercury/nf_scrna_qc:v3"
    }
    publishDir  path: "${params.outdir}/doublets",
                mode: "${params.copy_mode}",
                overwrite: "true"

    input:
        tuple(
            val(experiment_id),
            path(all_doublet_files)
        )

    output:
        // path("plots/*.pdf") optional true
        path("*.png") optional true
        // path("${experiment_id}__DoubletDecon_doublets_singlets.tsv"), emit: results
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
        
    main:
        // Identify multiplets using scrublet.
        gex_h5ad = CONVERT_MTX_TO_H5AD(channel__file_paths_10x).gex_h5ad
        input_channel = Channel.of()
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
        input_channel2.subscribe { println "1:: input_channel.out.citeseq_rsd: $it" }
        
        MERGE_DOUBLET_RESULTS(input_channel2) 


    emit:
        scrublet_paths = MERGE_DOUBLET_RESULTS.out.result
}

