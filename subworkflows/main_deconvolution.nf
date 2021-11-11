nextflow.enable.dsl=2

// main deconvolution modules, common to all input modes:

include { CELLSNP } from "$projectDir/modules/nf-core/modules/cellsnp/main"
include { SUBSET_GENOTYPE } from '../modules/nf-core/modules/subset_genotype/main'
include { VIREO } from '../modules/nf-core/modules/vireo/main'
include { GUZIP_VCF } from '../modules/nf-core/modules/guzip_vcf/main'
include { SOUPORCELL } from '../modules/nf-core/modules/souporcell/main'
include { SPLIT_DONOR_H5AD } from '../modules/nf-core/modules/split_donor_h5ad/main'
include { PLOT_DONOR_CELLS } from '../modules/nf-core/modules/plot_donor_cells/main'

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

                }else{
                    log.info "---We are using a full genotype input for Vireo----"
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

            if (params.vireo.run){


                SPLIT_DONOR_H5AD(vireo_out_sample_donor_ids.combine(ch_experiment_filth5, by: 0))

                // collect file paths to h5ad files in tsv tables:
                SPLIT_DONOR_H5AD.out.donors_h5ad_tsv
                .collectFile(name: "donors_h5ad.tsv",
                        newLine: false, sort: true,
                        seed: "experiment_id\tdonor\th5ad_filepath\n",
                        storeDir:params.outdir)

                // paste experiment_id and donor ID columns with __ separator
                SPLIT_DONOR_H5AD.out.exp__donors_h5ad_tsv
                .collectFile(name: "exp__donors_h5ad.tsv",
                        newLine: false, sort: true,
                        seed: "experiment_id\th5ad_filepath\n",
                        storeDir:params.outdir)

                SPLIT_DONOR_H5AD.out.donors_h5ad_assigned_tsv
                .collectFile(name: "donors_h5ad_assigned.tsv",
                        newLine: false, sort: true,
                        seed: "experiment_id\tdonor\th5ad_filepath\n",
                        storeDir:params.outdir)

                // paste experiment_id and donor ID columns with __ separator
                out_h5ad =SPLIT_DONOR_H5AD.out.exp__donors_h5ad_assigned_tsv
                .collectFile(name: "exp__donors_h5ad_assigned.tsv",
                        newLine: false, sort: true,
                        seed: "experiment_id\th5ad_filepath\n",
                        storeDir:params.outdir)
                
                SPLIT_DONOR_H5AD.out.h5ad_tsv
                .collectFile(name: "cellranger_as_h5ad.tsv",
                        newLine: true, sort: true, // only one line in each file to collate, without ending new line character, so add it here.
                        seed: "experiment_id\th5ad_filepath", // don't need \n here since newLine: true
                        storeDir:params.outdir)

                // all vireo() outputs collected -> plot_donor_ncells():
                vireo_out_sample_summary_tsv
                .collectFile(name: "vireo_donor_n_cells.tsv",
                        newLine: false, sort: true,
                        seed: "experiment_id\tdonor\tn_cells\n",
                        storeDir:params.outdir)
                .set{ch_vireo_donor_n_cells_tsv} // donor column: donor0, .., donorx, doublet, unassigned

                // paste experiment_id and donor ID columns with __ separator
                vireo_out_sample__exp_summary_tsv = vireo_out_sample__exp_summary_tsv
                .collectFile(name: "vireo_exp__donor_n_cells.tsv",
                        newLine: false, sort: true,
                        seed: "experiment_id\tn_cells\n",
                        storeDir:params.outdir)

                PLOT_DONOR_CELLS(ch_vireo_donor_n_cells_tsv)
                

        }else{
            out_h5ad ='None'
        }
        
    emit:
        out_h5ad
        vireo_out_sample__exp_summary_tsv

}
