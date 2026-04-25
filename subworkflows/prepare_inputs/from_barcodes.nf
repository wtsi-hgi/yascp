nextflow.enable.dsl=2

include {prep_collectmetadata; merge_metadata} from "$projectDir/modules/nf-core/modules/merge_metadata/main"

workflow from_barcodes {
    // Here we prepeare ann the required channels for the Cellbender and Deconvolution pipelines. 
    // QC steps and Data Handover bite on the merged h5ad created by MERGE_SAMPLES process and results directory output folder respectively.
    take: 
        channel_input_data_table
    main:
        log.info "... Prepearing inputs based on the 10x folder required for downstream analysis..."
        channel_input_data_table
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.experiment_id, row.n_pooled)}
            .set{ch_experiment_npooled}

        channel_input_data_table
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.experiment_id, "${row.data_path_10x_format}/possorted_genome_bam.bam" ,row.data_path_10x_format+'/filtered_feature_bc_matrix/barcodes.tsv.gz')}
            .set{pre_ch_experiment_bam_barcodes}

        channel__file_paths_10x =  channel_input_data_table
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row -> tuple(
            row.experiment_id,
            file("${row.data_path_10x_format}/filtered_feature_bc_matrix/barcodes.tsv.gz"),
            file("${row.data_path_10x_format}/filtered_feature_bc_matrix/features.tsv.gz"),
            file("${row.data_path_10x_format}/filtered_feature_bc_matrix/matrix.mtx.gz")
        )}

        channel__metadata =  channel_input_data_table
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row -> tuple(
            row.experiment_id,
            file("${row.data_path_10x_format}/metrics_summary.csv")
        )}
        
        prep_collectmetadata(channel__metadata)
        channel__metadata=merge_metadata(prep_collectmetadata.out.metadata.collect())


        channel_input_data_table
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row -> tuple(row.experiment_id, file(row.data_path_10x_format+'/raw_feature_bc_matrix'))}.set{ch_experimentid_paths10x_raw}

        channel_input_data_table
            .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row -> tuple(row.experiment_id, file(row.data_path_10x_format+'/filtered_feature_bc_matrix'))}.set{ch_experimentid_paths10x_filtered}


        if (params.split_h5ad_per_donor.run) {
            log.info "params.split_h5ad_per_donor.run=true: create pre_ch_experiment_filth5 from params.input_tables_column_delimiter"
            channel_input_data_table
                .splitCsv(header: true, sep: params.input_tables_column_delimiter)
            .map{row->tuple(row.experiment_id, "${row.data_path_10x_format}/filtered_feature_bc_matrix")}
            .set{pre_ch_experiment_filth5} // this channel is used for task 'split_donor_h5ad'
        }
        else {
            // create dummy channel
            log.info "using dummy channel for pre_ch_experiment_filth5, because unused."
            pre_ch_experiment_filth5 = Channel.from("foo").map { foo -> tuple("foo1","foo2") }
        }

        if (params.run_with_genotype_input) {
            log.info "You have selected params.vireo.run_with_genotype_input=true -> will run Vireo with genotype input. Input VCF and list of donors per experiment_id gathered from params.input_n_pooled_table)"
            if (params.genotype_input.subset_genotypes){
                log.info("----We will subset genotypes to the donors listed in the donor_vcf_ids for use in the Deconvolution----")
                // If we subset the genotypes, th ids in the donor_vcf_ids will be used to generate an individual vcf per pool.
                channel_input_data_table
                    .splitCsv(header: true, sep: params.input_tables_column_delimiter)
                    .map{row->tuple(row.experiment_id, params.genotype_input.full_vcf_file,params.genotype_input.full_vcf_file+'.csi', row.donor_vcf_ids)}
                    .set{pre_ch_experiment_donorsvcf_donorslist}
            }else{
                log.info('----We are using full VCF for each donor without subsetting----')
                // When we do not subset the genotypes, a full vcf will be fed in Vireo and soupocell assignments, but the number of pooled individuals will be used to detect the correct number of individuals in pool.
                // in this case we have to be aware that the last number is a number of donors pooled instead of IDs as per above
                channel_input_data_table
                    .splitCsv(header: true, sep: params.input_tables_column_delimiter)
                    .map{row->tuple(row.experiment_id, params.genotype_input.full_vcf_file,params.genotype_input.full_vcf_file+'.csi', row.n_pooled)}
                    .set{pre_ch_experiment_donorsvcf_donorslist}
            }

        } else {
            // Create dummy channel since we are not using VCFs, hence no genotype input
            log.info "We are not using genotypes in this pipeline run"
            pre_ch_experiment_donorsvcf_donorslist = Channel.from("foo").map { foo -> tuple("foo1","foo2","foo3",'foo4') }
        }

        pre_ch_experiment_filth5
            .set{ch_experiment_filth5}

        pre_ch_experiment_bam_barcodes
            .map { a,b,c -> tuple(a, file(b), file("${b}.bai"), file(c))}
            .set {ch_experiment_bam_bai_barcodes}

        pre_ch_experiment_donorsvcf_donorslist
            .set{ch_experiment_donorsvcf_donorslist}
        // ch_experiment_donorsvcf_donorslist.view()

    emit:
        ch_experiment_bam_bai_barcodes
        ch_experiment_npooled
        ch_experiment_filth5
        ch_experiment_donorsvcf_donorslist
        ch_experimentid_paths10x_raw
        ch_experimentid_paths10x_filtered
        channel__file_paths_10x
        channel__metadata
}
