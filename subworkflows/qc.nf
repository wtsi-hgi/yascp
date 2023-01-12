
// Load base.config by default for all pipelines - typically included in the nextflow config.
// Modules to include.
include {OUTLIER_FILTER} from "$projectDir/modules/nf-core/modules/outlier_filter/main"
include {PLOT_STATS} from "$projectDir/modules/nf-core/modules/plot_stats/main"
include {ESTIMATE_PCA_ELBOW} from "$projectDir/modules/nf-core/modules/estimate_pca_elbow/main"
include {SUBSET_PCS} from "$projectDir/modules/nf-core/modules/subset_pcs/main"
include {NORMALISE_AND_PCA} from "$projectDir/modules/nf-core/modules/normalise_and_pca/main"
include {HARMONY} from "$projectDir/modules/nf-core/modules/harmony/main"
include {BBKNN} from "$projectDir/modules/nf-core/modules/bbknn/main"
include {ADD_EXTRA_METADATA_TO_H5AD} from "$projectDir/modules/nf-core/modules/adata_manipulations/main"
include {LISI} from "$projectDir/modules/nf-core/modules/lisi/main"
include {UMAP; UMAP as UMAP_HARMONY; UMAP as UMAP_BBKNN;} from "$projectDir/modules/nf-core/modules/umap/main"
include {CLUSTERING; CLUSTERING as CLUSTERING_HARMONY; CLUSTERING as CLUSTERING_BBKNN;} from "$projectDir/modules/nf-core/modules/clustering/main"
include {CELL_HARD_FILTERS} from "$projectDir/modules/nf-core/modules/cell_hard_filters/main"


params.output_dir = "nf-qc_cluster"
    
workflow qc {
    take:
        file__anndata_merged
        file__cells_filtered
        assignments_all_pools
    main:
        log.info "--- Running QC metrics --- "
        // if(params.extra_metadata!=''){
        //     log.info '''--- Adding extra metadata to h5ad---'''
        //     ADD_EXTRA_METADATA_TO_H5AD(file__anndata_merged,params.extra_metadata)
        //     file__anndata_merged = ADD_EXTRA_METADATA_TO_H5AD.out.file__anndata
        // }else{
        //     log.info '''--- No extra metadata to add to h5ad ---'''
        // }
        file__anndata_merged.map{val1 -> tuple('full', val1)}.set{out1}

        CELL_HARD_FILTERS(file__anndata_merged,params.hard_filters_file,params.hard_filters_drop)
        if(params.hard_filters_file != "no_file__file_sample_qc"){
            file__anndata_merged = CELL_HARD_FILTERS.out.anndata
        }
        
        // Next we define an input channel to outlier filtering strategy in case if params.skip_preprocessing.gt_match_based_adaptive_qc_exclusion_pattern !=''
        // i.e - if we want to exclude a particular cohort that has been matched by gt match from the adaptive qc we feed this in the outlier_filter()
        if(params.skip_preprocessing.gt_match_based_adaptive_qc_exclusion_pattern !=''){
            gt_outlier_input = assignments_all_pools
        }else{
            gt_outlier_input = Channel.from("$projectDir/assets/fake_file.fq")
        }


        //FILTERING OUTLIER CELLS
        if (params.sample_qc.cell_filters.filter_outliers.run_process) {
            log.info """---Running automatic outlier cell filtering.----"""
            OUTLIER_FILTER(
                params.output_dir,
                file__anndata_merged,
                file__cells_filtered,
                params.sample_qc.cell_filters.filter_outliers.metadata_columns,
                params.sample_qc.cell_filters.filter_outliers.method,
                params.sample_qc.cell_filters.filter_outliers.outliers_fraction,
                params.sample_qc.cell_filters.filter_outliers.max_samples,
                params.anndata_compression_opts,
                gt_outlier_input,
                params.skip_preprocessing.gt_match_based_adaptive_qc_exclusion_pattern
            )
            file__anndata_merged = OUTLIER_FILTER.out.anndata
            file__cells_filtered = OUTLIER_FILTER.out.cells_filtered
        }

        

        NORMALISE_AND_PCA(params.output_dir+'/clustering',
            file__anndata_merged,
            params.mode,
            params.layer,
            params.genes_exclude_hvg,
            params.genes_score,
            params.reduced_dims.vars_to_regress.value)

        ESTIMATE_PCA_ELBOW(
            NORMALISE_AND_PCA.out.outdir,
            NORMALISE_AND_PCA.out.anndata,
            params.reduced_dims.n_dims.add_n_to_estimate
        )
        LI = ''

        if (params.reduced_dims.n_dims.auto_estimate) {
            log.info "n_pcs = automatically estimated."
            n_pcs = ESTIMATE_PCA_ELBOW.out.auto_elbow
        } else {
            log.info "n_pcs = Channel.from(params.reduced_dims.n_dims.value)"
            n_pcs = Channel.from(params.reduced_dims.n_dims.value)
        }

        SUBSET_PCS(
            NORMALISE_AND_PCA.out.outdir,
            NORMALISE_AND_PCA.out.anndata,
            NORMALISE_AND_PCA.out.metadata,
            NORMALISE_AND_PCA.out.pcs,
            NORMALISE_AND_PCA.out.param_details,
            n_pcs
        )

        PLOT_STATS(file__anndata_merged,
                    file__cells_filtered,
                    SUBSET_PCS.out.outdir,
                    SUBSET_PCS.out.anndata,
                    n_pcs)

        file__anndata_merged = NORMALISE_AND_PCA.out.anndata
        
        LI4 = PLOT_STATS.out.LI

        if (params.cluster.known_markers.run_process) {
            channel__cluster__known_markers = Channel
                .fromList(params.cluster.known_markers.value)
                .map{row -> tuple(row.file_id, file(row.file))}
        } else {
            channel__cluster__known_markers = tuple('', '')
        }

        // "Correct" PCs using Harmony or BBKNN
        if (params.harmony.run_process) {
            HARMONY(
                NORMALISE_AND_PCA.out.outdir,
                NORMALISE_AND_PCA.out.anndata,
                NORMALISE_AND_PCA.out.metadata,
                NORMALISE_AND_PCA.out.pcs,
                NORMALISE_AND_PCA.out.param_details,
                n_pcs,
                Channel.fromList( params.harmony.variables_and_thetas.value)
            )

            UMAP_HARMONY(
                HARMONY.out.outdir,
                HARMONY.out.anndata,
                HARMONY.out.metadata,
                HARMONY.out.pcs,
                HARMONY.out.reduced_dims,
                "False",
                params.umap.n_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value,
                params.umap.colors_quantitative.value,
                params.umap.colors_categorical.value,
                'harmony'
            )

            cluster_harmony__outdir = UMAP_HARMONY.out.outdir
            cluster_harmony__anndata = UMAP_HARMONY.out.anndata
            anndata =  UMAP_HARMONY.out.anndata
            outdir = UMAP_HARMONY.out.outdir
            cluster_harmony__metadata = UMAP_HARMONY.out.metadata
            cluster_harmony__pcs = UMAP_HARMONY.out.pcs
            cluster_harmony__reduced_dims = UMAP_HARMONY.out.reduced_dims
            

            CLUSTERING_HARMONY(
                cluster_harmony__outdir,
                cluster_harmony__anndata,
                cluster_harmony__metadata,
                cluster_harmony__pcs,
                cluster_harmony__reduced_dims,
                "False",  // use_pcs_as_reduced_dims
                params.cluster.number_neighbors.value,
                params.cluster.methods.value,
                params.cluster.resolutions.value,
                params.cluster.variables_boxplot.value,
                channel__cluster__known_markers,
                params.cluster_validate_resolution.sparsity.value,
                params.cluster_validate_resolution.train_size_cells.value,
                params.cluster_marker.methods.value,
                params.umap.n_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value,
                params.sccaf.min_accuracy         
            )
            
            lisi_input2 = HARMONY.out.reduced_dims_params.collect()
            LI1 = CLUSTERING_HARMONY.out.dummy_output
                
        }else{
            lisi_input2 = Channel.of([1, 'dummy'])
            LI1 = Channel.of([1, 'dummy'])
        }

        if (params.bbknn.run_process) {
            BBKNN(
                NORMALISE_AND_PCA.out.outdir,
                NORMALISE_AND_PCA.out.anndata,
                NORMALISE_AND_PCA.out.metadata,
                NORMALISE_AND_PCA.out.pcs,
                NORMALISE_AND_PCA.out.param_details,
                n_pcs,
                params.bbknn.batch_variable.value
            )

            UMAP_BBKNN(
                BBKNN.out.outdir,
                BBKNN.out.anndata,
                BBKNN.out.metadata,
                BBKNN.out.pcs,
                BBKNN.out.reduced_dims,
                "True",  // Don't look at the reduced_dims parameter
                ["-1"],  // params.cluster.number_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value,
                params.umap.colors_quantitative.value,
                params.umap.colors_categorical.value,
                'bbknn'
            )

            cluster_bbknn__outdir = UMAP_BBKNN.out.outdir
            cluster_bbknn__anndata = UMAP_BBKNN.out.anndata
            outdir = UMAP_BBKNN.out.outdir
            anndata = UMAP_BBKNN.out.anndata
            cluster_bbknn__metadata = UMAP_BBKNN.out.metadata
            cluster_bbknn__pcs = UMAP_BBKNN.out.pcs
            cluster_bbknn__reduced_dims = UMAP_BBKNN.out.reduced_dims

            CLUSTERING_BBKNN(
                cluster_bbknn__outdir,
                cluster_bbknn__anndata,
                cluster_bbknn__metadata,
                cluster_bbknn__pcs,
                cluster_bbknn__reduced_dims,
                "True",  // use_pcs_as_reduced_dims
                ["-1"],  // params.cluster.number_neighbors.value,
                params.cluster.methods.value,
                params.cluster.resolutions.value,
                params.cluster.variables_boxplot.value,
                channel__cluster__known_markers,
                params.cluster_validate_resolution.sparsity.value,
                params.cluster_validate_resolution.train_size_cells.value,
                params.cluster_marker.methods.value,
                ["-1"],  // params.umap.n_neighbors.value,
                params.umap.umap_init.value,
                params.umap.umap_min_dist.value,
                params.umap.umap_spread.value,
                params.sccaf.min_accuracy
            )
            
            lisi_input3 = BBKNN.out.reduced_dims_params.collect()
            LI2 = CLUSTERING_BBKNN.out.dummy_output
                
        }else{
            lisi_input3 = Channel.of([1, 'dummy'])
            LI2 = Channel.of([1, 'dummy'])
        }

        if (params.lisi.run_process) {
            lisi_input = SUBSET_PCS.out.reduced_dims_params.collect()
            lisi_input_first = lisi_input.mix(lisi_input2)
            lisi_input_second = lisi_input_first.mix(lisi_input3)

            LISI(
                NORMALISE_AND_PCA.out.outdir,
                NORMALISE_AND_PCA.out.metadata,
                params.lisi.variables.value,
                lisi_input_second.collect()
            )
            
            LI3 = LISI.out.outdir
        }else{
            LI3 = Channel.of([1, 'dummy'])
        }
        LI=LI1.mix(LI2)
        LI=LI.mix(LI3)
        LI=LI.mix(LI4)
    emit:
        LI
        file__anndata_merged
        
}
