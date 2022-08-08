
include {
    cluster;plot_phenotype_across_clusters;serialize_known_markers;plot_known_markers; cluster_validate_resolution_keras; plot_resolution_validate; cluster_markers; cellex_cluster_markers; prep_cellxgene;cluster_validate_resolution_sklearn
} from "./functions.nf"

include {umap_calculate_and_plot} from "../umap/functions.nf"

include {SCCAF} from '../sccaf/main'

workflow CLUSTERING {
    take:
        outdir
        anndata
        metadata
        pcs
        reduced_dims
        use_pcs_as_reduced_dims
        cluster__number_neighbors
        cluster__methods
        cluster__resolutions
        cluster__boxplot_variables
        cluster__known_markers
        cluster_validate_resolution__sparsity
        cluster_validate_resolution__train_size_cells
        cluster_marker__methods
        n_neighbors
        umap_init
        umap_min_dist
        umap_spread
        sccaf_minacc
    main:
        // Cluster the results, varying the resolution ------------------------
        cluster(
            outdir,
            anndata,
            metadata,
            pcs,
            reduced_dims,
            cluster__number_neighbors,
            cluster__methods,
            cluster__resolutions
        )

        // Boxplot of phenotype across clusters -------------------------------
        plot_phenotype_across_clusters(
            cluster.out.outdir,
            cluster.out.anndata,
            cluster__boxplot_variables
        )

        // // Dotplot of marker genes across clusters ----------------------------
        // // Serialize marker genes
        serialize_known_markers(
            cluster__known_markers
        )
        plot_known_markers(
            cluster.out.outdir,
            cluster.out.anndata,
            serialize_known_markers.out.marker_file
        )

        // // Validate the resolution
        // // Do not use cluster_validate_resolution_sklearn process.
        // cluster_validate_resolution_sklearn(
        //     cluster.out.outdir,
        //     cluster.out.anndata,
        //     cluster.out.metadata,
        //     cluster.out.pcs,
        //     cluster.out.reduced_dims,
        //     cluster.out.clusters,
        //     cluster_validate_resolution__sparsity,
        //     cluster_validate_resolution__train_size_cells
        // )
        if (params.utilise_gpu){
            if (params.cluster_validate_resolution_keras){

            
            cluster_validate_resolution_keras( 
                cluster.out.outdir,
                cluster.out.anndata,
                cluster.out.metadata,
                cluster.out.pcs,
                cluster.out.reduced_dims,
                cluster.out.clusters,
                cluster_validate_resolution__sparsity,
                cluster_validate_resolution__train_size_cells,
                cluster.out.outdir__reduced_dims
            )

            plot_resolution_validate(
                cluster_validate_resolution_keras.out.plot_input.groupTuple()
            )
            }
            
        }

        SCCAF(cluster.out.outdir,
          cluster.out.anndata,
          cluster.out.clusters,
          sccaf_minacc)




        // // Generate UMAPs of the results.
        umap_calculate_and_plot(
            cluster.out.outdir,
            cluster.out.anndata,
            cluster.out.pcs,
            cluster.out.reduced_dims,
            use_pcs_as_reduced_dims,
            "",
            "cluster",
            n_neighbors,
            umap_init,
            umap_min_dist,
            umap_spread
        )
        dummy_output=umap_calculate_and_plot.out.dummy_output
        // // Find marker genes for clusters
        cluster_markers(
            cluster.out.outdir,
            cluster.out.anndata,
            cluster.out.metadata,
            cluster.out.pcs,
            cluster.out.reduced_dims,
            cluster.out.clusters,
            cluster_marker__methods
        )

        // // Find marker genes for clusters using CELLEX
        cellex_cluster_markers(
            cluster.out.outdir,
            cluster.out.anndata
        )

        // Prep adata file for cellxgene website
        prep_cellxgene(
            cluster.out.outdir,
            cluster.out.anndata
        )

        emit:
            dummy_output

}