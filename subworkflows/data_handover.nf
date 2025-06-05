include { GATHER_DATA;  SPLIT_DATA_BY_STUDY} from "$projectDir/modules/local/gather_data/main"
include { ENCRYPT_DIR; ENCRYPT_TARGET } from "$projectDir/modules/local/encrypt"
include { TRANSFER;SUMMARY_STATISTICS_PLOTS } from "$projectDir/modules/local/summary_statistics_plots/main"
include { split_bam_by_donor } from "$projectDir/modules/local/cellranger_bam_per_donor"
// include { SUBSET_BAM_PER_BARCODES } from "$projectDir/modules/local/subset_bam_per_barcodes_and_variants/main"

workflow data_handover{
    take:
        outdir
        input_channel
        qc_input
        ch_poolid_csv_donor_assignments
        sample_possorted_bam_vireo_donor_ids
        genome
    main:
        log.info 'running data handover'

        if (params.gather_and_calculate_stats){
          GATHER_DATA(outdir,qc_input.collect(),input_channel)
          gh_out  = GATHER_DATA.out.outfiles_dataset
        }else{
          // Sometimes we do not want to calculate the statistics in depth
          gh_out  = Channel.from("$projectDir/assets/fake_file2.fq")
        }

        if (params.split_bam){
            GATHER_DATA.out.barcodes_files.flatten().map{sample -> tuple("${sample}".replaceFirst(/.*\//,"").replaceFirst(/\..*/,""),"${sample}".replaceFirst(/.*\//,"").replaceFirst(/\.tsv.*/,""),sample)}.set{barcodes}
            barcodes.combine(sample_possorted_bam_vireo_donor_ids, by: 0).set{full_split_chanel_input}
        } else {
          cram_encrypted_dirs = Channel.empty()
        }

        if (params.encrypt){
            ENCRYPT_DIR(GATHER_DATA.out.outfiles_dataset)
        
            cram_encrypted_dirs
              .collect()
              .set { ch_cram_encrypted_dirs }

            SPLIT_DATA_BY_STUDY(
              outdir,
              ENCRYPT_DIR.out.encrypted_dir,
              ch_cram_encrypted_dirs,
              ch_poolid_csv_donor_assignment_gathered
              )
        }

        SUMMARY_STATISTICS_PLOTS(outdir,gh_out,params.input_data_table)
        TRANSFER(SUMMARY_STATISTICS_PLOTS.out.summary_plots,params.rsync_to_web_file,outdir)

}
