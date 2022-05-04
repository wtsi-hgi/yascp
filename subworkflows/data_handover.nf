include { GATHER_DATA } from '../modules/nf-core/modules/gather_data/main'
include { SUMMARY_STATISTICS_PLOTS } from '../modules/nf-core/modules/summary_statistics_plots/main'
include { split_bam_by_donor } from "../modules/local/cellranger_bam_per_donor"

workflow data_handover{
    take:
        outdir
        qc_input
        
    main:
        log.info 'running data handover'
        GATHER_DATA(outdir,qc_input.collect(),params.input_data_table)
        SUMMARY_STATISTICS_PLOTS(outdir,GATHER_DATA.out.outfiles_dataset,params.input_data_table)
        if (params.split_bam){
            split_bam_by_donor(main_deconvolution.out.sample_possorted_bam_vireo_donor_ids)
        }
        
        // We should also generate a report simmilar to the MultiQC that can be easily shared around in a html format.
        // If we run it in sanger through the gilab CI we transfer the data to the local website for data investigations.
        
}
