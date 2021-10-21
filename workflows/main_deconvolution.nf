nextflow.enable.dsl=2

// main deconvolution modules, common to all input modes:

include { CELLSNP } from "$projectDir/modules/nf-core/modules/cellsnp/main"
include { SUBSET_GENOTYPE } from '../modules/nf-core/modules/subset_genotype/main'
include { VIREO } from '../modules/nf-core/modules/vireo/main'
include { GUZIP_VCF } from '../modules/nf-core/modules/guzip_vcf/main'
include { SOUPORCELL } from '../modules/nf-core/modules/souporcell/main'

workflow  main_deconvolution {

    take:
		ch_experiment_bam_bai_barcodes
		ch_experiment_npooled
		ch_experiment_filth5
		ch_experiment_donorsvcf_donorslist

    main:
		log.info "#### running workflow main_deconvolution() ..."


	if (params.run_with_genotype_input) {
		if (params.genotype_input.subset_genotypes){
			log.info "---We are subsetting genotypes----"

			SUBSET_GENOTYPE(ch_experiment_donorsvcf_donorslist.map { experiment, donorsvcf, donorslist -> tuple(experiment,
							file(params.cellsnp.vcf_candidate_snps),
							file(donorsvcf),
							donorslist)})
		}

	}

	ch_experiment_donorsvcf_donorslist.map { experiment, donorsvcf, donorslist -> tuple(experiment, donorslist.replaceAll(/,/, " ").replaceAll(/"/, ""))}.set{donors_in_lane}


    // cellsnp() outputs -> vireo():
	if (params.vireo.run){

        // cellsnp() from pipeline provided inputs:
        CELLSNP(ch_experiment_bam_bai_barcodes,
            Channel.fromPath(params.cellsnp.vcf_candidate_snps).collect())

		// Vireo:
		if (params.run_with_genotype_input) {
			log.info "---running Vireo with genotype input----"
			// for each experiment_id to deconvolute, subset donors vcf to its donors and subset genomic regions.
			if (params.genotype_input.subset_genotypes){
				log.info "---We are using subset genotypes running Vireo----"
                CELLSNP.out.cellsnp_output_dir.combine(ch_experiment_npooled, by: 0)
                    .combine(SUBSET_GENOTYPE.out.samplename_subsetvcf, by: 0).set{full_vcf}
                full_vcf.view()
			}else{
				log.info "---We are using a full genotype input for Vireo----"

				// cellsnp.out.cellsnp_output_dir.combine(cellsnp.out.cellsnp_output_dir, by:0).view()
                CELLSNP.out.cellsnp_output_dir.combine(ch_experiment_npooled, by: 0).set{full_vcf}
				full_vcf.map { experiment, cellsnpvcf, npooled -> tuple(experiment,cellsnpvcf,npooled,file(params.genotype_input.full_vcf_file))}.set {full_vcf}

			}


		}
		// Vireo without genotype input:
		else {
			log.info "-----running Vireo without genotype input----"
            CELLSNP.out.cellsnp_output_dir.combine(ch_experiment_npooled, by: 0).set{full_vcf}
            full_vcf.map {experiment, cellsnp, npooled -> tuple(experiment, cellsnp, npooled,[])}.set{full_vcf}

		}

        VIREO(full_vcf)
        vireo_out_sample_summary_tsv = VIREO.out.sample_summary_tsv
        vireo_out_sample__exp_summary_tsv = VIREO.out.sample__exp_summary_tsv
        vireo_out_sample_donor_ids = VIREO.out.sample_donor_ids
	}


	if (params.souporcell.run){
		// This runs the Souporcell Preprocessing
		if (params.souporcell.use_raw_barcodes) {
			// read raw cellranger barcodes per pool for souporcell
			channel.fromPath(params.souporcell.path_raw_barcodes_table)
					.splitCsv(header: true, sep: params.input_tables_column_delimiter)
				.map{row->tuple(row.experiment_id, row.data_path_barcodes.replaceFirst(/${params.replace_in_path_from}/, params.replace_in_path_to))}
				.set{ch_experiment_rawbarcodes}
			ch_experiment_bam_bai_barcodes
				.map {a,b,c,d -> tuple(a,b,c)}
				.combine(ch_experiment_rawbarcodes, by: 0)
				.combine(ch_experiment_npooled, by: 0)
				.set {ch_experiment_bam_bai_barcodes_npooled}
		} else if (! params.souporcell.use_raw_barcodes) {
			ch_experiment_bam_bai_barcodes
				.combine(ch_experiment_npooled, by: 0)
				.set{ch_experiment_bam_bai_barcodes_npooled}
		}

		// This runs the Souporcell
		if (params.run_with_genotype_input) {

			if (params.genotype_input.subset_genotypes){
				log.info "---We are using subset genotypes running Suporcell----"
				// here combine the ch_experiment_bam_bai_barcodes_npooled with the output of subset
				GUZIP_VCF(SUBSET_GENOTYPE.out.samplename_subsetvcf)
				// ch_experiment_bam_bai_barcodes_npooled.combine(SUBSET_GENOTYPE.out.samplename_subsetvcf, by: 0).set{full_vcf}
				ch_experiment_bam_bai_barcodes_npooled.combine(GUZIP_VCF.out.souporcell_vcf, by: 0).set{full_vcf}


			}else{
				log.info "---We are using a full genotype input for Suporcell----"
				// this however currently doesnt work and the individuals have to be provided.

				// here just add the full vcf path to each of the ch_experiment_bam_bai_barcodes_npooled
				GUZIP_VCF(tuple('full_vcf', file(params.genotype_input.full_vcf_file)))
				GUZIP_VCF.out.souporcell_vcf.map { sample, vcf -> vcf }.set{vcf_file}
				// we would now flatten this and take the the file of the tuple to be used next.
				ch_experiment_bam_bai_barcodes_npooled.map {
					samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters ->
					tuple(samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters)
					}.set{full_vcf}
				full_vcf.combine(vcf_file).set{full_vcf}

			}
			full_vcf.combine(donors_in_lane, by: 0).set{full_vcf}
		}
		else{
			log.info "-----running Suporcell without genotype input----"
			// here make add an empty entry [] to the ch_experiment_bam_bai_barcodes_npooled
			ch_experiment_bam_bai_barcodes_npooled.map {
				samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters ->
				tuple(samplename, bam_file, bai_file, barcodes_tsv_gz, souporcell_n_clusters,[],[])
				}.set{full_vcf}
		}


		// Now that channel is created run suporcell
		SOUPORCELL(full_vcf,
			Channel.fromPath(params.souporcell.reference_fasta).collect())

	}


}
