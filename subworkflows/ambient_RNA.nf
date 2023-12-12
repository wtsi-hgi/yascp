
// Load base.config by default for all pipelines - typically included in the nextflow config.
include { CELLBENDER } from '../modules/nf-core/modules/cellbender/main'
include { SPLIT_CITESEQ_GEX; DSB } from '../modules/nf-core/modules/citeseq/main'
include {capture_cellbender_files} from "$projectDir/modules/nf-core/modules/cellbender/functions"

workflow ambient_RNA {
    take:
        ch_experimentid_paths10x_raw
		ch_experimentid_paths10x_filtered
        channel__metadata
        ab_data
    main:
        log.info params.input_data_table
        log.info """---Running Cellbender pipeline ---"""

        capture_cellbender_files(params.cellbender_location,"${params.outdir}/nf-preprocessing",params.input_data_table)
        capture_cellbender_files.out.alt_input.flatten().map{sample -> tuple("${sample}".replaceFirst(/.*\/captured\//,"").replaceFirst(/\/.*/,""),sample)}.set{alt_input}
        capture_cellbender_files.out.cb_to_use_downstream.flatten().map{sample -> tuple("${sample}".replaceFirst(/.*\/cellbender\//,"").replaceFirst(/\/.*/,""),sample)}.set{cb_Filtered_pre2}
        ch_experimentid_paths10x_raw.join(alt_input, by: [0], remainder: true).set{post_ch_experimentid_paths10x_raw}
        ch_experimentid_paths10x_filtered.join(alt_input, by: [0], remainder: true).set{post_ch_experimentid_paths10x_filtered}
        post_ch_experimentid_paths10x_raw.filter{ it[1] != null }.filter{ it[2] == null }.map{row -> tuple(row[0], row[1])}.set{ch_experimentid_paths10x_raw_2}
        post_ch_experimentid_paths10x_filtered.filter{ it[1] != null }.filter{ it[2] == null }.map{row -> tuple(row[0], row[1])}.set{ch_experimentid_paths10x_filtered_2}
        post_ch_experimentid_paths10x_raw.filter{ it[1] != null }.filter{ it[2] != null }.map{row -> tuple(row[0],row[2])}.set{alt_input3}

        CELLBENDER(ch_experimentid_paths10x_raw_2,ch_experimentid_paths10x_filtered_2,channel__metadata)
        
        cb_Filtered_pre = CELLBENDER.out.cellbender_downstream
        // If we rerun this with a provided cb path this will not emit anything. 
        cb_Filtered =cb_Filtered_pre2.concat(cb_Filtered_pre)
        // AMBIENTNESS_QUANTIFICATION()
        // cb_Filtered_pre2 = capture_cellbender_files.out.cellbender_downstream
        // cb_Filtered.subscribe { println "cb_Filtered: $it" }
        if (params.citeseq){
            // Here we capture the filtered matrix and applyt the DSB normalisation on the citeseq datasets
            cb_Filtered.join(ab_data, by: [0], remainder: false).set{cb_ab}
            cb_ab.join(ch_experimentid_paths10x_raw, by: [0], remainder: false).set{cb_ab_raw}
            // DSB(cb_ab_raw)
        }

        alt_input1 = CELLBENDER.out.cellbender_path
        cellbender_path = alt_input3.concat(alt_input1)
        // cellbender_path=CELLBENDER.out.cellbender_path
    emit:
        // results_list
        cellbender_path
}
