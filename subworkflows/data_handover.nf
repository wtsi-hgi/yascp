include { GATHER_DATA } from '../modules/nf-core/modules/gather_data/main'

workflow data_handover{
    take:
        outdir
        cellbender_input
        qc_input
        
        
    main:
        log.info 'running data handover'
        GATHER_DATA(outdir,cellbender_input)
        
        SUMMARY_STATISTICS_PLOTS(outdir,cellbender_input,GATHER_DATA.out.outfiles_dataset2)
        // We also generate a report.

        // If we run it in sanger we transfer the data to the local website.
        
}