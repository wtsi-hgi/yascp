params {
     /////////////////////////
        // QC inputs:
     /////////////////////////


    file_paths_10x = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/mo11_work/outputs/test_qc/inputs/file_paths_10x.tsv"

    mode = "conventional_h5ad"
    // tsv table with columns "experiment_id" "donor" and "h5ad_filepath"
    file_paths_h5ad = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/mo11_work/outputs/test_qc/inputs/file_paths_h5ad.tsv"

    // Tab-delimited file containing sample metadata.
    file_metadata = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/mo11_work/outputs/test_qc/inputs/file_metadata.tsv"

    //    --genes_score Tab-delimited file with genes to use to score cells. Must contain ensembl_gene_id and score_id columns. If one score_id == "cell_cycle", then requires a grouping_id column with "G2/M" and "S". If no filter, then pass an empty file.
    genes_score = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/mo11_work/outputs/test_qc/inputs/genes_score_v001.tsv"

    //    --genes_exclude_hvg Tab-delimited file with genes to exclude from highly variable gene list. Must contain ensembl_gene_id column. If no filter, then pass an empty file.
    genes_exclude_hvg = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/mo11_work/outputs/test_qc/inputs/genes_remove_hvg_v001.tsv"

    //    --output_dir Directory name to save results to. (Defaults to 'nf-qc_cluster')
    output_dir = "${projectDir}/../results_nf_qc_cluster"

    run_multiplet = true

    // data_handover{
    //     run_process = true
    //     output_dir = 'Franke_with_genotypes_data_handover'
    //     //cellranger_raw_files_table = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/fetch/Connect_Val_2/results/raw.Submission_Data_Pilot_UKB.file_paths_10x.tsv"
    //     cellranger_raw_files_table = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/deconv/franke_data_postcellbender/inputs/franke_data.Pilot_3.raw.file_paths_10x_farm_paths.tsv"
    //     cellbender_files_table = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/cellbender/franke_data/outputs/qc_cluster_input_files/file_paths_10x-cellbender_params__epochs_250__learnrt_1pt0Eneg7__zdim_100__zlayer_500__lowcount_10-FPR_0pt05.tsv"
    //     deconvolution_dir = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/deconv/Franke_with_genotypes/results/split_donor_h5ad"
    //     Fetch_path = ''
    //     Cellbender_path = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/cellbender/franke_data'
    //     Deconvolution_path = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/deconv/Franke_with_genotypes'
    //     deconvolution_files_table =  "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/deconv/Franke_with_genotypes/results/exp__donors_h5ad.tsv"

    // }


    // transfer_to_web{
    //     run_process = true
    //     //# destination: 'ubuntu@172.27.22.139:/volume/scRNA_test_app/scrna_static_and_media_files/bin'
    //     destination = 'ubuntu@172.27.22.139:/volume/scRNA_test_app/scrna_static_and_media_files/media'
    //     tranche{
    //         tranche_run_name = 'Franke_with_genotypes' //its important that the name is unique, if not this will overwrite the existing run
    //         Fetch_path = ''
    //         Cellbender_path = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/cellbender/franke_data'
    //         Deconvolution_path = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/deconv/Franke_with_genotypes'
    //     }
    // }


    /////////////////////////
    // QC input finish
    /////////////////////////


    input = 'cellbender'
    cellsnp_input_table_mode = "from_barcodes"//"from_barcodes"
    // input_data_table colums: experiment_id   data_path_10x_format
    // hgi note: from nf_cellbender pipeline
    utilise_gpu = true
    ///////////////////////
    // Cellbender output
    input_data_table = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/mo11_work/outputs/test_qc/inputs/vcf_names_pooled.tsv'

    input_tables_column_delimiter = '\t' // set 'tsv' or 'csv': whether input_data_table (and other input tables) have tab-separted ('tsv') or comma-separated columns ('csv').

    celltypist { // cf. https://github.com/Teichlab/celltypist
        run = false // whether to run 'celltypist' task

        // run celltypist on filtered barcodes straight from cellranger outputs:
        //   must have columns 'experiment_id' (ID of cellranger run) and 'data_path_filt_h5d' (absolute path to cellranger filtered filtered_feature_bc_matrix.h5)
        // path_filtered_h5 = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/fetch/wbc_mult_donor/results/Submission_Data_Pilot_UKB.file_paths_10x.tsv'

        remove_workdir = false // // whether to remove all work dirs of this task when workflow{} is finished.
        copy_mode = "rellink" // choose "rellink", "symlink", "move" or "copy".
        // Make sure copy_mode is either "copy" or "move" when remove_workdir = true

        // specify models to use,
        //   from default available models from 'models.download_models(force_update = True)'
        //   cf. https://github.com/Teichlab/celltypist
        //   comma separated list of celltypist default models to use:
        models = ['Immune_All_High.pkl','Immune_All_Low.pkl',
            'Immune_Blood_High.pkl','Immune_Blood_Low.pkl']
    }

    run_with_genotype_input=true
	genotype_input {
        subset_genotypes = true
	    // path to donor vcfs:
        full_vcf_file = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/deconv/Franke_with_genotypes/inputs/franke_gt.vcf.gz'
    }

    cellsnp {
        run = true // whether to run 'cellsnp' task
        remove_workdir = false // // whether to remove all work dirs of this task when workflow{} is finished.
        copy_mode = "rellink" // choose "rellink", "symlink", "move" or "copy".
        // Make sure copy_mode is either "copy" or "move" when remove_workdir = true

        vcf_candidate_snps = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/nf_ci_deconv_inputs/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz"
        // this list of candidate SNPs for cellSNP comes from link at https://github.com/single-cell-genetics/cellSNP
        // i.e., https://sourceforge.net/projects/cellsnp/files/SNPlist/genome1K.phase3.SNP_AF5e2.chr1toX.hg38.vcf.gz/download

        // cellSNP CLI parameters:
        min_maf = "0.1" // cellSNP --minMAF
        min_count = "60" // cellSNP --minCOUNT
        p = "20" // cellSNP -p
    }

    vireo {
        // choose either run=true or run_with_genotype_input=true (to run Vireo without or with genotype input VCF).
        run = true // whether to run 'vireo' task
        remove_workdir = false // // whether to remove all work dirs of this task when workflow{} is finished.
        copy_mode = "rellink" // choose "rellink", "symlink", "move" or "copy".
        // Make sure copy_mode is either "copy" or "move" when remove_workdir = true

        // with_genotype_input = true // set to true or false. If true, will need the following two inputs:
    }

    plot_donor_ncells {
        run = false // whether to run 'plot_donor_ncells' task (multi-pages pdf for all samples Vireo deconvolutions)
        remove_workdir = false // // whether to remove all work dirs of this task when workflow{} is finished.
        copy_mode = "rellink" // choose "rellink", "symlink", "move" or "copy".
        // Make sure copy_mode is either "copy" or "move" when remove_workdir = true
        plotnine_dpi = "100"
    }

    split_h5ad_per_donor {
        run = true // whether to run 'split_h5ad_donor' task
        remove_workdir = false // // whether to remove all work dirs of this task when workflow{} is finished.
        copy_mode = "rellink" // choose "rellink", "symlink", "move" or "copy".
        // Make sure copy_mode is either "copy" or "move" when remove_workdir = true

        // next, optional arguments for python script split_h5ad_per_donor.py (cf. bin directory)
        input_h5_genome_version = "GRCh38"
        print_modules_version = "True"
        plot_n_cells_per_vireo_donor = "True"
        write_donor_level_filtered_cells_h5 = "True"
        plotnine_dpi = "100"
        anndata_compression_level = "6" // with gzip: compression_opts sets the compression level and may be an integer from 0 to 9, default is 4.

        // temporary fix:
        // specify manually absolute /lustre path to nextflow results directory
        // this will be used for generating absolute paths to h5ad deconvoluted objects,as listed in output tsv tables..
        // must append / with \\ as this string will be used by sed replacement.
        absolute_results_path = "\\/lustre\\/scratch123\\/hgi\\/projects\\/ukbb_scrna\\/pipelines\\/Pilot_UKB\\/deconv\\/Franke_with_genotypes\\/results"
        // absolute_results_path = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/deconv/Franke_with_genotypes/results"
    }

    subset_genotype {
        remove_workdir = true // // whether to remove all work dirs of this task when workflow{} is finished.
        copy_mode = "rellink" // choose "rellink", "symlink", "move" or "copy".
        // Make sure copy_mode is either "copy" or "move" when remove_workdir = true
    }
    souporcell {
        run = true // whether to run 'souporcell' task

        // run souporcell on raw or filtered barcodes:
        use_raw_barcodes = false

        // provides RAW cellranger barcodes as input to souporcell:
        //  must have columns experiment_id and data_path_barcodes
        path_raw_barcodes_table = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/fetch/Franke_with_genotypes/results/raw.Submission_Data_Pilot_UKB.file_paths_10x.tsv'

        remove_workdir = false // // whether to remove all work dirs of this task when workflow{} is finished.
        copy_mode = "rellink" // choose "rellink", "symlink", "move" or "copy".
        // Make sure copy_mode is either "copy" or "move" when remove_workdir = true

        // to get this ref file, cf. script inputs/get_souporcell_ref.sh
        reference_fasta = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/pipeline_inputs/deconv/refdata-cellranger-GRCh38-3.0.0/fasta/genome.fa'
        // path is secure lustre path from openstack instance, not Sanger farm path
        //   reference_fasta = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/pipeline_inputs/deconv/refdata-cellranger-GRCh38-3.0.0'
    }

    plot_souporcell_vs_vireo {
        // if both souporcell and vireo are set to run,
        //   then a Venn diagram can be made to compare cell assignements (to deconvoluted donor/cluster, doublet or unassigned)
        run = false // whether to run 'venn_diagram_souporcell_vs_vireo' task
        remove_workdir = false // // whether to remove all work dirs of this task when workflow{} is finished.
        copy_mode = "rellink" // choose "rellink", "symlink", "move" or "copy".
    }

    // other input parameters common to all input modes:

    // The following is used if any module above has replace_in_path = true.
    // use if runnning in Sanger FCE/Openstack (not necessary for Sanger LSF farm) :
    replace_in_path = false // if true, will replace the /lustre path to data_path_10x_format for FCE mapping of /lustre:
    replace_in_path_from = '/lustre/scratch123/hgi/mdt1/projects/ukbb_scrna' // part of path in real /lustre
    replace_in_path_to = '/lustre/scratch123' // to replace for Openstack mapping of /lustre path


    // TODO: the following two options are not fully coded yet, so setting them to true won't do anything (as of Jan 26th 2021).
    // the following are for one-off tasks run after workflow completion to clean-up work directories:
    on_complete_remove_workdirs = false // will remove work dirs (effectively un-caching) of selected tasks even if completed successfully. Make sure that copy_mode is also set to copy or move.
    on_complete_remove_workdir_failed_tasks = false // will remove work dirs of failed tasks (.exitcode file not 0)

}
