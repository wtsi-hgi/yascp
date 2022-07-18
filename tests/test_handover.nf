#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {data_handover} from "${projectDir}/subworkflows/data_handover"

params.output_dir = "results_test"
params.input_data_table = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/hp3_dev/elgh_yascp/inputs.tsv"

workflow {
  main:
    data_handover()
    generate_emitting_channel(ch_azimuth)
    generate_emitting_channel.out.oufn
      .collect()
      .set{ch_gathering_trigger}
    ch_gathering_trigger.subscribe onNext: { println "gathering", it }, onComplete: { println 'Done gathering' }

//    test_channel(
//      ch_gathering_trigger,
//      params.output_dir
//      )
    handover("${workDir}/../${params.output_dir}",qc_finish_dummy)
  }
