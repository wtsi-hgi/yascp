params {

    lisi{
        run_process=true
    }
    replace_genotype_ids=false
    write_h5=true
    cluster_validate_resolution_keras = true
    // run_celltype_assignment = true
    project_name = 'T_Cell_Bio_Response'
    filter_outliers = false
    extra_sample_metadata =""
    output_dir = outdir= "${launchDir}/recluster_resolutions"
    cellex_cluster_markers=true 
    cluster_markers = false
    normalise_andata = false
    skip_handover = true
    // output_dir = outdir= "${launchDir}/results"
    // run_celltype_assignment=true
    split_ad_per_bach=true //if not splitting the celltype assignment will be run on full tranche
    // input_data_table = "$outdir/handover/Summary_plots/$RUN_ID/Fetch Pipeline/Input/input_table.tsv"
    // cellbender_location="${output_dir}/nf-preprocessing/cellbender" //!!!!! if cellbender is run already then can skip this by selecting  input = 'existing_cellbender' instead input = 'cellbender'
    // existing_cellsnp="${output_dir}/cellsnp"
    cellbender_location="/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/harriet/qc/results_11_09_2023/nf-preprocessing/cellbender" //!!!!! if cellbender is run already then can skip this by selecting  input = 'existing_cellbender' instead input = 'cellbender'
    existing_cellsnp="/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/harriet/qc/results/cellsnp"

    skip_preprocessing = true
    // file__anndata_merged = '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/harriet_analysis/230313_hb58_yascp_analysis/231114_h5ad_files_for_MCC/231120_TCs_only_regressed_counts_HVGs.h5ad'

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
    // hard_filters_file       = "${projectDir}/../sample_qc.yml"
    // hard_filters_drop = false //#This indicates whether we want to drop the cells that fail hard filters of just flag them

    cluster{
        description = """Parameters for clustering. All pairwise combinations of
        method and resolution will be performed."""
        number_neighbors{
            description = """Number of neighbors. If <= 0, uses number of unique
            experiment_id."""
            value = 15
        }
        methods{
            description = 'Clustering method. Valid options [leiden|louvain].'
            value = 'leiden'
        }
        resolutions{
            description = 'Clustering resolution.'
            value = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        }

        variables_boxplot{
            decription = 'Generate boxplots of these variables for each cluster.'
            value ='n_cells,total_counts,pct_counts_gene_group__mito_transcript'
        }

        known_markers{
            run_process = false
            description = """Files with markers that will be used to generate
            dotplots. Each marker file should be the full path and have the
            following columns: cell_type, hgnc_symbol. The following columns
            are optional: p_value_adj. Use "" for a single entry in the
            file_id and file value to indicate no plots."""
            value = [
                [ file_id: 'SmillieCS_31348891', file: '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/marker_gene_db-raw_data/database/celltypes/colon/SmillieCS-31348891/database.tsv' ],
                [ file_id: 'ParikhK_30814735', file: '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/marker_gene_db-raw_data/database/celltypes/colon/ParikhK-30814735/database.tsv' ],
                [ file_id: 'JamesKR_32066951', file: '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/marker_gene_db-raw_data/database/celltypes/colon-immune/JamesKR-32066951/database.tsv' ]
            ]
        }




    }  
    bbknn{
        run_process = true
    }

    celltype_assignment{
        run_celltype_assignment=false
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
        memory = 300.GB
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
