include { GATHER_DATA } from '../modules/nf-core/modules/gather_data/main'
include { SUMMARY_STATISTICS_PLOTS } from '../modules/nf-core/modules/summary_statistics_plots/main'

workflow data_handover{
    take:
        outdir
        qc_input
        
    main:
        log.info 'running data handover'
        GATHER_DATA(outdir,qc_input,params.input_data_table)
        SUMMARY_STATISTICS_PLOTS(outdir,GATHER_DATA.out.outfiles_dataset,params.input_data_table)
        // We should also generate a report simmilar to the MultiQC that can be easily shared around in a html format.
        // If we run it in sanger we transfer the data to the local website for data investigations.
        
}
