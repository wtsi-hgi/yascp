include { GATHER_DATA;  SPLIT_DATA_BY_STUDY} from '../modules/nf-core/modules/gather_data/main'
include { ENCRYPT_DIR } from '../modules/local/encrypt'
include { TRANSFER;SUMMARY_STATISTICS_PLOTS } from '../modules/nf-core/modules/summary_statistics_plots/main'
include { split_bam_by_donor } from "../modules/local/cellranger_bam_per_donor"

workflow data_handover{
    take:
        outdir
        qc_input
        ch_poolid_csv_donor_assignments


    main:
        log.info 'running data handover'
        GATHER_DATA(outdir,qc_input.collect(),params.input_data_table)
        SPLIT_DATA_BY_STUDY(outdir, GATHER_DATA.out.outfiles_dataset, ch_poolid_csv_donor_assignments)
        if (params.encrypt){
            ENCRYPT_DIR(SPLIT_DATA_BY_STUDY.out.outdir_ukbb)
        }

        if (params.split_bam){
            split_bam_by_donor(main_deconvolution.out.sample_possorted_bam_vireo_donor_ids)
        }

        SUMMARY_STATISTICS_PLOTS(outdir,GATHER_DATA.out.outfiles_dataset,params.input_data_table)
        // We also generate a report.
        // If we run it in sanger we transfer the data to the local website.
        TRANSFER(SUMMARY_STATISTICS_PLOTS.out.summary_plots)

}
