params {
    RUN = "$RUN_ID"
    PATH2 = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/fetch/${RUN}"
    input = 'existing_cellbender' 
    lisi{
        run_process=false
    }
    replace_genotype_ids=false
    webtransfer = true
    write_h5=true
    encrypt = false
    // hard_filters_file = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/Blood_Fresh/yascp/conf/extra_confs/sample_qc.yml" //# This may point to the sample_qc.yml input which will apply hard filters to the merged cells.
    rsync_to_web_file='/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/fetch/SETUP_fech_input_prep/scripts/freeze1/rsync_to_web.sh'
    cohorts_to_drop_from_GT_Relatednes_check='GT_UKBB'
    project_name = 'Cardinal_Freeze1'
    cellbender_location='/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/ELGH_GT_All_Rerun/results_rsync2/results/nf-preprocessing/cellbender' //!!!!! if cellbender is run already then can skip this by selecting  input = 'existing_cellbender' instead input = 'cellbender'
    cellbender_resolution_to_use='0pt1'
    extra_metadata = "$PATH2/results/yascp_inputs/Extra_Metadata.tsv"
    extra_sample_metadata ="$PATH2/results/yascp_inputs/Extra_Metadata_Donors.tsv"
    genotype_phenotype_mapping_file = '/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/secret/bridge.txt'
    output_dir = outdir= "${launchDir}/results_ct_adaptiveqc/results"
        
    // output_dir = outdir= "${launchDir}/results"
    run_celltype_assignment=true
    split_ad_per_bach=true //if not splitting the celltype assignment will be run on full tranche
    input_data_table = "$PATH2/results/yascp_inputs/input.tsv" //this has to be a full path
    // input_data_table = "$outdir/handover/Summary_plots/$RUN_ID/Fetch Pipeline/Input/input_table.tsv"
    //!!!!! if cellbender is run already then can skip this by selecting  input = 'existing_cellbender' instead input = 'cellbender'
    existing_cellsnp="/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/qc/Blood_Fresh/results_hard_filters/cellsnp"
    existing_vireo="$output_dir/deconvolution/infered_genotypes"
    reference_assembly_fasta_dir_bam_split = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/genotypes/10x_reference_assembly_prefix"
    // output_dir = outdir= "${launchDir}/results_vireo_inp/results"
    webtransfer = true
	genotype_input {
        run_with_genotype_input=true
        vireo_with_gt=false
        posterior_assignment = false //if this is set to true, we will perform the genotype donor matching after the deconvolution is performed.
        subset_genotypes = false //if activated this will use th IDs provided to estimate the donors, otherwise it will match against full cohort
        full_vcf_file = '' //this could be a list of vcfs, in which case have to merge them 
        tsv_donor_panel_vcfs = "/lustre/scratch123/hgi/projects/ukbb_scrna/pipelines/Pilot_UKB/genotypes/donor_panels_vcf_paths.tsv"
    }
    souporcell {
        run = false
    }
    skip_preprocessing{
        value=true
        gt_match_file="$output_dir/deconvolution/deconvolution_results/vireo_gt_fix/assignments_all_pools.tsv" //We prvide this if we want to exclude a particular samples matched to a ceirtain GT cohortc from the adaptive qc
        gt_match_based_adaptive_qc_exclusion_pattern = 'U937;THP1' //We run the adaptive QC on these patterns independently regardless on assigned celltype.
        file__anndata_merged = "$output_dir/merged_h5ad/1.pre_QC_adata.h5ad"
        file__cells_filtered = "$output_dir/merged_h5ad/pre_QC_adata-cell_filtered_per_experiment.tsv.gz"

    }
    harmony{
        run_process= true
    }
    umap{
        run_process = true
        colors_quantitative{
            description = 'Comma separated string of quantitative variables that will be used to color points.'
            value = 'n_cells,total_counts,pct_counts_gene_group__mito_transcript,prob_doublet,pct_counts_gene_group__ribo_rna,Azimuth:predicted.celltype.l2.score,Azimuth:mapping.score,log10_ngenes_by_count'
        }
        colors_categorical{
            description = 'Comma separated string of categorical variables that will be used to color points.'
            value = 'cell_passes_qc,cell_passes_qc-per:Azimuth:L0_predicted.celltype.l2,experiment_id,Azimuth:predicted.celltype.l2,Celltypist:Immune_All_Low:predicted_labels,Celltypist:Immune_All_High:predicted_labels,donor_id'
        }
    }

    mads_categories ='pct_counts_gene_group__mito_transcript,pct_counts_gene_group__mito_protein,pct_counts_gene_group__ribo_protein,pct_counts_gene_group__ribo_rna,total_counts,n_genes_by_counts,log10_ngenes_by_count'
    hard_filters_file       = "${projectDir}/../sample_qc.yml"
    hard_filters_drop = false //#This indicates whether we want to drop the cells that fail hard filters of just flag them
       
    bbknn{
        run_process = true
    }

    celltype_assignment{
        run_celltype_assignment=true
        run_azimuth=true
        run_keras=false
        run_celltypist=true
    }
    reduced_dims{
        vars_to_regress{
            value = ''   
        }
    }

}

process {

    withName: plot_distributions{
        containerOptions = "--containall --cleanenv --workdir /tmp -B /tmp"
    }

    withName: cellex_cluster_markers{
        maxForks=7
        memory = 300.GB
    }
    
    withName: GATHER_DATA{
        maxForks=7
        memory = 100.GB
    }
    withName: LISI{
        maxForks=7
        memory = 500.GB
    }
    withName: cluster_validate_resolution_keras{
        memory = 300.GB
    }

    withName: umap_calculate_and_plot{
        memory = 300.GB
    }

    withName: sccaf_assess_clustering{
        memory = 300.GB
    }
    
}

singularity {
  enabled = true
  cacheDir   = "${baseDir}/singularity"
  runOptions = '--bind /lustre --no-home'
}