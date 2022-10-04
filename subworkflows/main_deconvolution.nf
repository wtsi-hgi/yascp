nextflow.enable.dsl=2

// main deconvolution modules, common to all input modes:

include { CELLSNP } from "../modules/nf-core/modules/cellsnp/main"
include { SUBSET_GENOTYPE } from '../modules/nf-core/modules/subset_genotype/main'
include { VIREO } from '../modules/nf-core/modules/vireo/main'
include { GUZIP_VCF } from '../modules/nf-core/modules/guzip_vcf/main'
include { SOUPORCELL } from '../modules/nf-core/modules/souporcell/main'
include { SPLIT_DONOR_H5AD } from '../modules/nf-core/modules/split_donor_h5ad/main'
include { PLOT_DONOR_CELLS } from '../modules/nf-core/modules/plot_donor_cells/main'
include {MULTIPLET} from "../modules/nf-core/modules/multiplet/main"
include {SOUPORCELL_VS_VIREO} from "../modules/nf-core/modules/plot_souporcell_vs_vireo/main"
include { REPLACE_GT_ASSIGNMENTS_WITH_PHENOTYPE; ENHANCE_VIREO_METADATA_WITH_DONOR } from '../modules/nf-core/modules/genotypes/main'
include { match_genotypes } from './match_genotypes'
include {REPLACE_GT_DONOR_ID } from '../modules/nf-core/modules/genotypes/main'
workflow  main_deconvolution {

    take:
		ch_experiment_bam_bai_barcodes
		ch_experiment_npooled
		ch_experiment_filth5
		ch_experiment_donorsvcf_donorslist
        channel__file_paths_10x

    main:
		log.info "#### running DECONVOLUTION workflow #####"
        if (params.run_with_genotype_input) {
            if (params.genotype_input.posterior_assignment){
                log.info('Deconvolution will run without genotypes, but Genotypes will be used after deconvolution to assign correct labels')
            }else{
                if (params.genotype_input.subset_genotypes){
                    log.info "---We are subsetting genotypes----"

                    SUBSET_GENOTYPE(ch_experiment_donorsvcf_donorslist.map { experiment, donorsvcf,donortbi, donorslist -> tuple(experiment,
                                    file(donorsvcf),file(donortbi),
                                    donorslist)})
                }
            }


        }
        ch_poolid_donor_assignment = Channel.empty()
        ch_experiment_donorsvcf_donorslist.map { experiment, donorsvcf, donorstbi,donorslist -> tuple(experiment, donorslist.replaceAll(/,/, " ").replaceAll(/"/, ""))}.set{donors_in_lane}

        CELLSNP(ch_experiment_bam_bai_barcodes,
            Channel.fromPath(params.cellsnp.vcf_candidate_snps).collect())

        MULTIPLET(
            params.output_dir,
            channel__file_paths_10x,
            params.sample_qc.cell_filters.filter_multiplets.expected_multiplet_rate,
            params.sample_qc.cell_filters.filter_multiplets.n_simulated_multiplet,
            params.sample_qc.cell_filters.filter_multiplets.multiplet_threshold_method,
            params.sample_qc.cell_filters.filter_multiplets.scale_log10
        )

        // cellsnp() outputs -> vireo():
        if (params.vireo.run){
            // Here we run Vireo software to perform the donor deconvolution. Note that we have coded the pipeline to be capable in using
            // the full genotypes as an input and also subset to the individuals provided as an input in the donor_vcf_ids column. The
            // VIREO:
            if (params.run_with_genotype_input) {
                log.info "---running Vireo with genotype input----"
                // for each experiment_id to deconvolute, subset donors vcf to its donors and subset genomic regions.
                if (params.genotype_input.subset_genotypes){
                    // Here we subset the genotypes. This happens if the input.nf contains subset_genotypes = true
                    log.info "---We are using subset genotypes running Vireo----"
                    CELLSNP.out.cellsnp_output_dir.combine(ch_experiment_npooled, by: 0)
                        .combine(SUBSET_GENOTYPE.out.samplename_subsetvcf, by: 0).set{full_vcf}
                }else{
                    // Here we do not subset the genotypes and match against the full cohort provided as an input. This happens if subset_genotypes = false
                    log.info "---We are using a full genotype input for Vireo----"
                    CELLSNP.out.cellsnp_output_dir.combine(ch_experiment_npooled, by: 0).set{full_vcf}
                    full_vcf.map { experiment, cellsnpvcf, npooled -> tuple(experiment,cellsnpvcf,npooled,file(params.genotype_input.full_vcf_file),file(params.genotype_input.full_vcf_file+'.csi'))}.set {full_vcf}
                }
            }
            // Vireo without genotype input:
            else {
                // Vireo can also be run without the genotypes, and the performance is equally good then running with.
                // here we run it without the genotypes and the donors are labeled as donor0, donor1 etc, dependant on the number of donors set in the input file n_pooled column
                log.info "-----running Vireo without genotype input----"
                CELLSNP.out.cellsnp_output_dir.combine(ch_experiment_npooled, by: 0).set{full_vcf}
                full_vcf.map {experiment, cellsnp, npooled -> tuple(experiment, cellsnp, npooled,[],[])}.set{full_vcf}
            }

            // When all the channels are prpeared accordingly we exacute the vireo with the prpeared channel.
            full_vcf.filter { experiment, cellsnp, npooled, t,ti -> npooled != '1' }.set{full_vcf2}
            full_vcf.filter { experiment, cellsnp, npooled, t,ti -> npooled == '1' }.set{not_deconvoluted}
            VIREO(full_vcf2)
            vireo_out_sample_donor_vcf = VIREO.out.infered_vcf
            vireo_out_sample_summary_tsv = VIREO.out.sample_summary_tsv
            vireo_out_sample__exp_summary_tsv = VIREO.out.sample__exp_summary_tsv
            vireo_out_sample_donor_ids = VIREO.out.sample_donor_ids
        }


        if (params.souporcell.run){
            // YASCP pipeline is also capable in running SOUPORCELL instead of VIREO. If activated SOUPORCELL will be used.
            // yascp currently doesnt have an option to take souporcell assignments as downstream instead of vireo but this will be added shortly.
            // This runs the Souporcell Preprocessing

            if (params.souporcell.use_raw_barcodes) {
                // read raw cellranger barcodes per pool for souporcell
                log.info """Here use the raw barcodes"""
            } else if (! params.souporcell.use_raw_barcodes) {
                ch_experiment_bam_bai_barcodes
                    .combine(ch_experiment_npooled, by: 0)
                    .set{ch_experiment_bam_bai_barcodes_npooled}
            }

            // This runs the Souporcell
            // Similarly to the VIREO Soupocell can be run with and without genotypes and folowing prpeares the inputs accordingly to each option.
            // Soupocell cant digest a gz file gence we extract the data.
            if (params.run_with_genotype_input) {

                if (params.genotype_input.subset_genotypes){
                    // this will run the soupocell with the subset genotypes. This happens if the input.nf contains subset_genotypes = true
                    log.info "---We are using subset genotypes running Suporcell----"
                    GUZIP_VCF(SUBSET_GENOTYPE.out.samplename_subsetvcf)
                    // ch_experiment_bam_bai_barcodes_npooled.combine(SUBSET_GENOTYPE.out.samplename_subsetvcf, by: 0).set{full_vcf}
                    ch_experiment_bam_bai_barcodes_npooled.combine(GUZIP_VCF.out.souporcell_vcf, by: 0).set{full_vcf}

                }else{
                    log.info "---We are using a full genotype input for Suporcell----"
                    // this however currently doesnt work and the individuals have to be provided.
                    // here just add the full vcf path to each of the ch_experiment_bam_bai_barcodes_npooled
                    GUZIP_VCF(tuple('full_vcf', file(params.genotype_input.full_vcf_file)))
                    GUZIP_VCF.out.souporcell_vcf.map { sample, vcf -> vcf }.set{vcf_file}
                    ch_experiment_bam_bai_barcodes_npooled.map {
                        samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters ->
                        tuple(samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters)
                        }.set{full_vcf}
                    full_vcf.combine(vcf_file).set{full_vcf}

                }
                full_vcf.combine(donors_in_lane, by: 0).set{full_vcf}
            }
            else{
                // Soupocell can also run without the genotypes. This will prpeare channels to run it withoutt.
                log.info "-----running Suporcell without genotype input----"
                // here make add an empty entry [] to the ch_experiment_bam_bai_barcodes_npooled
                ch_experiment_bam_bai_barcodes_npooled.map {
                    samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters ->
                    tuple(samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters,[],[])
                    }.set{full_vcf}
            }

            // When all the channels are prepeared then we can run the soupocell accordingly.
            SOUPORCELL(full_vcf,
                Channel.fromPath(params.souporcell.reference_fasta).collect())
            // Regardless if the Soupocell is run with or without genotypes we still need to match the donor ids with the cluster ids since this does not happen automatically is Soupocell.
        }

        if (params.run_with_genotype_input & params.genotype_input.posterior_assignment) {
            match_genotypes(vireo_out_sample_donor_vcf)

            out_gt = match_genotypes.out.donor_match_table
            ch_poolid_donor_assignment = match_genotypes.out.pool_id_donor_assignments_csv
            // //here we fix the genotype ids to the ones matched by the GT match similarly to what vireo would do.
            REPLACE_GT_DONOR_ID(VIREO.out.all_required_data , out_gt.collect())
            if (params.replace_genotype_ids){

                REPLACE_GT_DONOR_ID.out.sample_donor_vcf.set{vireo_out_sample_donor_vcf}
                REPLACE_GT_DONOR_ID.out.sample_summary_tsv.set{vireo_out_sample_summary_tsv}
                REPLACE_GT_DONOR_ID.out.sample__exp_summary_tsv.set{vireo_out_sample__exp_summary_tsv}
                REPLACE_GT_DONOR_ID.out.sample_donor_ids.set{vireo_out_sample_donor_ids}


            }
            REPLACE_GT_DONOR_ID.out.assignments
                    .collectFile(name: "assignments_all_pools.tsv",
                            newLine: false, sort: true,
                            keepHeader: true,
                            // skip:1,
                            storeDir:params.outdir+'/deconvolution/vireo_gt_fix')
            // otherwise we enhance the vireo metadata report with sample ids.

        }

        //here have to fix the vireo outputs based on the GT matching.
        if (params.replace_genotype_ids){
            REPLACE_GT_DONOR_ID(VIREO.out.all_required_data , ch_poolid_donor_assignment.collect())
            REPLACE_GT_DONOR_ID.out.sample_donor_vcf.set{vireo_out_sample_donor_vcf}
            REPLACE_GT_DONOR_ID.out.sample_summary_tsv.set{vireo_out_sample_summary_tsv}
            REPLACE_GT_DONOR_ID.out.sample__exp_summary_tsv.set{vireo_out_sample__exp_summary_tsv}
            REPLACE_GT_DONOR_ID.out.sample_donor_ids.set{vireo_out_sample_donor_ids}
        }


        if (params.vireo.run){

                // The folowing downstream tasks prepeares the plots, splits the donors according to vireo ids and generates the summary files.
                // These folowing outputs should be cleaned up and have been left as they were in the deconvolution pipeline initially.
                // We also need to account that others may want to run soupocell instead and this folowing is not capable in digesting the cluster ids by soupocell currently.

                ch_experiment_bam_bai_barcodes
                  .map { samplename, bam_file, bai_file, barcodes_tsv_gz -> tuple(samplename, file(bam_file)) }
                  .combine(vireo_out_sample_donor_ids, by: 0 )
                  .set { ch_experiment_bam_vireo_donor_ids }

                // If sample is not deconvoluted we will use scrublet to detect the doublets and remove them.
                not_deconvoluted.map{ experiment, donorsvcf, npooled,t -> tuple(experiment, 'None')}.set{not_deconvoluted2}
                file_cellmetadata = MULTIPLET.out.file__cellmetadata
                scrublet_paths = MULTIPLET.out.scrublet_paths

                // making 2 channels, 1 that is deconvoluted and another that isnt
                split_channel = vireo_out_sample_donor_ids.combine(ch_experiment_filth5, by: 0)
                split_channel2 = not_deconvoluted2.combine(ch_experiment_filth5, by: 0)

                // combining these 2 channels in one
                split_channel3 = split_channel.mix(split_channel2)

                // adding the scrublet paths to the channel.
                split_channel4 = split_channel3.combine(scrublet_paths, by: 0)

                split_channel5 = split_channel4.map{
                val_sample, val_donor_ids_tsv, val_filtered_matrix_h5, path_scrublet ->
                [  val_sample,
                file(val_donor_ids_tsv),
                file(val_filtered_matrix_h5),
                    path_scrublet,
                params.outdir]}

                SPLIT_DONOR_H5AD(split_channel5)

                // collect file paths to h5ad files in tsv tables:
                SPLIT_DONOR_H5AD.out.donors_h5ad_tsv
                .collectFile(name: "donors_h5ad.tsv",
                        newLine: false, sort: true,
                        seed: "experiment_id\tdonor\th5ad_filepath\n",
                        storeDir:params.outdir+'/deconvolution/filepaths')

                // paste experiment_id and donor ID columns with __ separator
                SPLIT_DONOR_H5AD.out.exp__donors_h5ad_tsv
                .collectFile(name: "exp__donors_h5ad.tsv",
                        newLine: false, sort: true,
                        seed: "experiment_id\th5ad_filepath\n",
                        storeDir:params.outdir+'/deconvolution/filepaths')

                SPLIT_DONOR_H5AD.out.donors_h5ad_assigned_tsv
                .collectFile(name: "donors_h5ad_assigned.tsv",
                        newLine: false, sort: true,
                        seed: "experiment_id\tdonor\th5ad_filepath\n",
                        storeDir:params.outdir+'/deconvolution/filepaths')

                // paste experiment_id and donor ID columns with __ separator
                out_h5ad =SPLIT_DONOR_H5AD.out.exp__donors_h5ad_assigned_tsv
                .collectFile(name: "exp__donors_h5ad_assigned.tsv",
                        newLine: false, sort: true,
                        seed: "experiment_id\th5ad_filepath\n",
                        storeDir:params.outdir+'/deconvolution/filepaths')

                SPLIT_DONOR_H5AD.out.h5ad_tsv
                .collectFile(name: "cellranger_as_h5ad.tsv",
                        newLine: true, sort: true, // only one line in each file to collate, without ending new line character, so add it here.
                        seed: "experiment_id\th5ad_filepath", // don't need \n here since newLine: true
                        storeDir:params.outdir+'/deconvolution/filepaths')

                // all vireo() outputs collected -> plot_donor_ncells():
                vireo_out_sample_summary_tsv
                .collectFile(name: "vireo_donor_n_cells.tsv",
                        newLine: false, sort: true,
                        seed: "experiment_id\tdonor\tn_cells\n",
                        storeDir:params.outdir+'/deconvolution/filepaths')
                .set{ch_vireo_donor_n_cells_tsv} // donor column: donor0, .., donorx, doublet, unassigned

                out_split = SPLIT_DONOR_H5AD.out.donor_n_cells
                // here before splitting we may want to add extra metadata if it is based on the genotype ids.



                // paste experiment_id and donor ID columns with __ separator
                vireo_out_sample__exp_summary_tsv = out_split
                .collectFile(name: "vireo_exp__donor_n_cells.tsv",
                        newLine: false, sort: true,
                        seed: "experiment_id\tn_cells\n",
                        storeDir:params.outdir+'/deconvolution/filepaths')


                if (params.run_with_genotype_input & params.genotype_input.posterior_assignment) {
                    if (!params.replace_genotype_ids & params.extra_sample_metadata!='' & params.add_donor_metadata){
                        // Here we have sample level metadata but we have chosen to keep the donor ids.
                        // in this scenario we enhance the donor level vireo metadata file to add the donor metadata to the h5ads and eventually to the donor and teanche report
                        ENHANCE_VIREO_METADATA_WITH_DONOR(params.extra_sample_metadata,vireo_out_sample__exp_summary_tsv,REPLACE_GT_DONOR_ID.out.assignments.collect())
                        vireo_out_sample__exp_summary_tsv = ENHANCE_VIREO_METADATA_WITH_DONOR.out.replaced_vireo_exp__donor_n_cells_out
                    }
                }


                PLOT_DONOR_CELLS(ch_vireo_donor_n_cells_tsv)


        } else{
            out_h5ad =Channel.fromPath(params.cellsnp.vcf_candidate_snps).collect()
            vireo_out_sample__exp_summary_tsv=Channel.fromPath(params.cellsnp.vcf_candidate_snps).collect()
            vireo_out_sample_donor_vcf = Channel.empty()
            ch_experiment_bam_vireo_donor_ids = Channel.empty()
        }

        if (params.souporcell.run && params.vireo.run) {
            // This is an orringinal comparison of Vireo and soupocell, this needs to be cleaned up.
            SOUPORCELL_VS_VIREO(
                vireo_out_sample_donor_ids
                .combine(SOUPORCELL.out.souporcell_output_files.map {a,b,c,d -> tuple(a,b)},
                    by: 0))
        }

    emit:
        out_h5ad
        vireo_out_sample__exp_summary_tsv
        vireo_out_sample_donor_vcf
        ch_poolid_csv_donor_assignments = ch_poolid_donor_assignment
        sample_possorted_bam_vireo_donor_ids = ch_experiment_bam_vireo_donor_ids
}
