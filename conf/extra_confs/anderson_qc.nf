params{

    help                 = false
    run_multiplet        = false
    mode                 = "conventional"
    layer                = "none"
    file_sample_qc       = "no_file__file_sample_qc"
    file_cellmetadata    = "no_file__file_cellmetadata"
    file_anndata         = "no_file__file_anndata"
    genes_exclude_hvg    = "no_file__genes_exclude_hvg"
    genes_score          = "no_file__genes_score"
    anndata_compression_opts = 9

    anndata_compression_opts = 9
    run_celltype_assignment = true
    genes_score = "${projectDir}/assets/genes_score_v001.tsv"
    genes_exclude_hvg = "${projectDir}/assets/genes_remove_hvg_v001.tsv"
    input = 'cellbender'
    metadata_key_column{
            description = """Column in metadata that matches the experiment_id column in
        tenx_data."""
            value = 'experiment_id'
    }

    celltype_prediction {
	    // # parameters for the keras celltype prediction python script
        keras {
	    keras_model = '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy/sc_qc_cluster/pca_plus5/cellbender_fpr0pt1_filteroutlier_0pt1/nf_results/normalize=total_count.vars_to_regress=none.hvg_exclude=data-variable_gene_filter.scores=data-gene_scores/reduced_dims-vars_to_regress=none-bbknn.batch=experiment_id.n_pcs=29/cluster.number_neighbors=-1.method=leiden.resolution=3.0/validate_resolution/adata-normalized_pca-bbknn-umap-clustered-sparsity_l1=1pt0E-4-train_size_cells=-1.h5'
	    keras_weights_df = '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/analysis/freeze_003/ti-cd_healthy/sc_qc_cluster/pca_plus5/cellbender_fpr0pt1_filteroutlier_0pt1/nf_results/normalize=total_count.vars_to_regress=none.hvg_exclude=data-variable_gene_filter.scores=data-gene_scores/reduced_dims-vars_to_regress=none-bbknn.batch=experiment_id.n_pcs=29/cluster.number_neighbors=-1.method=leiden.resolution=3.0/validate_resolution/adata-normalized_pca-bbknn-umap-clustered-sparsity_l1=1pt0E-4-train_size_cells=-1-weights.tsv.gz'
            h5_layer = 'log1p_cp10k'
	    // #cloned from github:
	    // # https://github.com/andersonlab/sc_nextflow-studies/blob/main/gut-freeze003/ti-cd_healthy/results/cluster_annotations/data-cluster_labels.csv
	    keras_model_cluster_labels = "${projectDir}/assets/keras_cluster/data-cluster_labels.csv"
	    filter_top_cell_probabilities = '0.5,0.75'
	    save_all_probabilities = '--save_all_probabilities'
	}
    }


    celltypist {
        run = true
        description = """https://github.com/Teichlab/celltypist"""
        remove_workdir = false
        copy_mode = "rellink"
        models = ['Immune_All_High.pkl','Immune_All_Low.pkl']
    }

    sample_qc{

        description = """Parameters for sample QC prior to merge.
            Filters are applied to all samples."""
        cell_filters{
            description = """Cell filters. Each bullet point is a seperate filter.
            Cells that evaluate to true for any of these filters, will be
            removed. Filters under 'all_samples' are applied to all samples.
            Filters under a sample id are applied to that specific sample."""
            filter_multiplets{
                description = """Parameters for scrublet. Runs prior to filters
                    below. Note scale_log10 should be 'True|False'.
                    Output from multiplet analysis will be added to the final
                    AnnData object. The flag only works if
                    file_cellmetadata from the main nextflow call is not set."""
                run_process= true
                expected_multiplet_rate= 0.1
                n_simulated_multiplet= 100000
                multiplet_threshold_method = 'threshold_li'
                scale_log10 = 'False'

            }
            all_samples{
                description = 'Cell filters applied to all samples.'
                value = 'pct_counts_gene_group__mito_transcript > 50 and n_genes_by_counts <= 100 and scrublet__predicted_multiplet == True'
            }



            filter_outliers{
                description= """After applying cell filters based on cutoffs, apply
                    an filter to remove outliers using cell information from
                    metadata_columns.
                    Recommended settings:
                        method: IsolationForest
                        metadata_columns: n_cells
                        outliers_fraction: 0.0
                        max_cells: 0.1 or 0.99
                    Notes:
                    * Filters applied across the full dataset (i.e., not on a per
                        sample level).
                    * Valid methods: LocalOutlierFactor, IsolationForest,
                        EllipticEnvelope, OneClassSVM.
                    * outliers_fraction: anticipated fraction of outlier cells.
                        If 0.0, then runs sklearn methods with 'auto' as the
                        anticipated number of outlier cells.
                    * max_samples: The fraction of cells to draw from X to train
                        each estimator. Only valid if method == IsolationForest."""
                run_process= true
                method= 'IsolationForest'
                metadata_columns= 'pct_counts_gene_group__mito_transcript,log1p_total_counts,log1p_n_genes_by_counts'
                outliers_fraction= 0.0
                max_samples= 0.1
            }




        }

        downsample_cells_fraction{
            description = """Downsample to this fraction of the number of
                observations (sc.pp.subsample). Use "" to indicate no downsampling.
                Example: if 0.8, will drop 20% of cells."""
            value= ''
        }

        downsample_cells_n{
            description = """Downsample to this number of observations
                (sc.pp.subsample). Use "" to indicate no downsampling. Example: if
                200 then 200 total cells will be kept."""
            value = ''
        }

        downsample_feature_counts{
            description= """Downsample the number of feature counts by this fraction.
                Use "" to indicate no downsampling."""
            value = ''
        }

    }

    plots_qc{

        description = """Parameters for basic QC plots after sample merge. All
        pairwise combinations of facet_columns and
        variable_columns_distribution_plots will be performed."""
        facet_columns{
            description= 'Column to facet all QC plots by.'
            value= 'experiment_id,sex,disease_status,bead_version,lane,ack_lysis_buffer'
        }

        variable_columns_distribution_plots{
            description= """Plot the distributions of these variables (histogram and
            ecdf)."""
            value = 'total_counts,pct_counts_gene_group__mito_transcript,pct_counts_gene_group__mito_protein,pct_counts_gene_group__ribo_protein'
        }

    }

    reduced_dims{
        run_downstream_analysis = false
        description = """Parameters for dimensionality reduction (principal component
        and harmony calculations). All pairwise combinations of
        vars_to_regress and n_dims will be performed."""

        vars_to_regress{
            description = """Comma separated string of variables to regress. Use "" to
            indicate no regession is to be performed."""
            //# value = 'pct_counts_gene_group__mito_transcript'
            value = ''
            
        }

        n_dims{
            description =  """Number of dimensions to use for calculating PCs,
            harmony corrected PCs, and clusters. Value should be int.
            If auto_estimate == true then the 'knee' of a scree plot is
            automatically estimated via the kneedle algorithm and used.
            If auto_estimate == true then value is ignored.
            add_n_to_estimate is added to the elbow estimate."""
            auto_estimate = true
            add_n_to_estimate = 5
            value = 27
        }

    }

    harmony{
        run_process= false
        description= 'Parameters for harmony.'
        variables_and_thetas{
            description = 'Tuples of metadata columns and corresponding thetas.'
            value = [
                [ variable: 'experiment_id', theta: 0.5 ],
                [ variable: 'experiment_id', theta: 1.0 ],
                [ variable: 'experiment_id', theta: 2.0 ],
                ]
        }
    }

    bbknn{
        run_process = true
        description = 'Parameters for BBKNN'
        batch_variable{
            description= 'Variable to use for batch correction'
            value = 'experiment_id'
        }
    }

    lisi{
        run_process= true
        description = 'Parameters for Local Inverse Simpsons Index (LISI).'
        variables{
            description = 'Metadata variables to compute LISI over.'
            value = 'experiment_id'
        }
    }

    cluster{
        description = """Parameters for clustering. All pairwise combinations of
        method and resolution will be performed."""
        number_neighbors{
            description = """Number of neighbors. If <= 0, uses number of unique
            experiment_id."""
            value = [30,50,100,150]
        }
        methods{
            description = 'Clustering method. Valid options [leiden|louvain].'
            value = 'leiden'
        }
        resolutions{
                        description = 'Clustering resolution.'
            value = [0.5,1.0,1.5,1.75,2.0,2.25,2.5,3.0,3.5,4.0,5.0]
        }

        variables_boxplot{
            decription = 'Generate boxplots of these variables for each cluster.'
            value ='total_counts,age,pct_counts_gene_group__mito_transcript,pct_counts_gene_group__mito_protein,pct_counts_gene_group__ribo_protein,scrublet__multiplet_zscores,S_score,G2M_score,b_cell__activation,b_cell__differentiation,colon_compartment__epithelial,colon_compartment__immune,colon_compartment__stromal,crypt__axis_score,t_cell__coinhibition_signature,t_cell__costimulation_signature,t_cell__cytotoxicity_signature,t_cell__cd8_cytotoxic,t_cell__cell_cycling,t_cell__central_memory,t_cell__immunoregulation,t_cell__resident_memory,villus__bottom_enterocytes,villus__top_enterocytes,hallmark_apoptosis,hallmark_hypoxia,mense_hypoxia_up,hallmark_inflammatory_response,inflammation_signature__smilliecs_31348891,colon_compartment__epithelial,colon_compartment__immune,colon_compartment__stromal'
        }

        known_markers{
            run_process = true
            description = """Files with markers that will be used to generate
            dotplots. Each marker file should be the full path and have the
            following columns: cell_type, hgnc_symbol. The following columns
            are optional: p_value_adj. Use "" for a single entry in the
            file_id and file value to indicate no plots."""
            value = [
                [ file_id: 'SmillieCS_31348891', file: '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/marker_gene_db-raw_data/database/celltypes/colon/SmillieCS-31348891/database.tsv' ],
                [ file_id: 'ParikhK_30814735', file: '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/marker_gene_db-raw_data/database/celltypes/colon/ParikhK-30814735/database.tsv' ],
                [ file_id: 'JamesKR_32066951', file: '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/marker_gene_db-raw_data/database/celltypes/colon-immune/JamesKR-32066951/database.tsv' ],
                [ file_id: 'KincheJ_30270042', file: '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/marker_gene_db-raw_data/database/celltypes/colon-mesenchyme/KincheJ-30270042/database.tsv' ],
                [ file_id: 'ElmentaiteR_20200206937110', file: '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/marker_gene_db-raw_data/database/celltypes/terminal_ileum-epithelial/ElmentaiteR-2020.02.06.937110/database.tsv' ],
                [ file_id: 'AndersonGroupVeli_201906', file: '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/marker_gene_db-raw_data/database/celltypes/terminal_ileum/AndersonGroup-Veli_2019_06/database.tsv' ],
                [ file_id: 'WangY_31753849', file: '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/marker_gene_db-raw_data/database/celltypes/colon-ileum-rectum/WangY-31753849/database.tsv' ],
                [ file_id: 'CorridoniD_32747828', file: '/lustre/scratch119/humgen/projects/sc-eqtl-ibd/data/marker_gene_db-raw_data/database/celltypes/colon-immune/CorridoniD-32747828/database.tsv' ],
            ]
        }




    }

    cluster_validate_resolution{
        description: 'Parameters for cluster resolution validation.'
        sparsity{
            description=  """LogisticRegression sparsity or inverse of regularization
            strength; must be a positive float. Like in support vector
            machines, smaller values specify stronger regularization."""
            value = 0.0001
        }
        train_size_cells{
            description= """Number of cells to use for training. Set to -1 for
                to use default of 2/3 of total cells. Do not change this parameter
                unless you know what you are doing."""
            value = -1
        }
    }

    cluster_marker{
        description = 'Parameters for identifying cluster marker features.'
        methods{
            description = 'Method for marker detection. Valid options [wilcoxon|logreg].'
                    value = 'wilcoxon'
        }
    }

    umap{
        run_process = true
        description = 'Parameters for umap.'
        colors_quantitative{
            description = 'Comma separated string of quantitative variables that will be used to color points.'
            value = 'total_counts,age,pct_counts_gene_group__mito_transcript,pct_counts_gene_group__mito_protein,pct_counts_gene_group__ribo_protein,scrublet__multiplet_zscores,S_score,G2M_score,b_cell__activation,b_cell__differentiation,colon_compartment__epithelial,colon_compartment__immune,colon_compartment__stromal,crypt__axis_score,t_cell__coinhibition_signature,t_cell__costimulation_signature,t_cell__cytotoxicity_signature,t_cell__cd8_cytotoxic,t_cell__cell_cycling,t_cell__central_memory,t_cell__immunoregulation,t_cell__resident_memory,villus__bottom_enterocytes,villus__top_enterocytes,hallmark_apoptosis,hallmark_hypoxia,mense_hypoxia_up,hallmark_inflammatory_response,inflammation_signature__smilliecs_31348891,colon_compartment__epithelial,colon_compartment__immune,colon_compartment__stromal'
        }
        colors_categorical{
            description = 'Comma separated string of categorical variables that will be used to color points.'
            value = 'experiment_id,sex,disease_status,scrublet__predicted_multiplet,bead_version,phase'
        }
        n_neighbors{
            description = """Number of neighbors for sc.pp.neighbors call.
                Recommended value between 2-100. If you expect each cell type
                cluster to be shared across all experiments/samples, then setting
                this number to the number of experiments/samples is a good place to
                start. Note: values separated with a comma will be run within the
                same script call (rather than swarm)."""
            value= [30,50,100]
        }

        umap_init{
            description = 'How to initialize the low dimensional embedding.'
            value = 'spectral'
        }

        umap_min_dist{
            description = """The effective minimum distance between embedded points.
                Recommended value between 0-1. Note: values separated with a comma
                will be run within the same script call (rather than swarm)."""
            value = [0.0,0.25,0.5,0.75,1.0]
        }

        umap_spread{
            description = """The effective scale of embedded points.
                Recommended value between 0-3. Note: values separated with a comma
                will be run within the same script call (rather than swarm)."""
            value = [0.5,1.0,2.5]
        }
    }

    sccaf{
        description = 'sccaf'
        run_assessment = true
        run_optimization = false
        min_accuracy = 0.92
        default_leiden_res = 4.0
    }

    azimuth{
        run_process = true
        celltype_refset = 'celltype.l2'
    }

}

