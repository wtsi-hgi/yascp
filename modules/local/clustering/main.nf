
include {
    CLUSTER;PLOT_PHENOTYPE_ACROSS_CLUSTERS;SERIALIZE_KNOWN_MARKERS;PLOT_KNOWN_MARKERS; CLUSTER_VALIDATE_RESOLUTION_KERAS; PLOT_RESOLUTION_VALIDATE; CLUSTER_MARKERS; CELLEX_CLUSTER_MARKERS; PREP_CELLXGENE
} from "./functions.nf"

include {UMAP_CALCULATE_AND_PLOT} from "../umap/functions.nf"

// include {SCCAF} from '../sccaf/main'

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
        Channel.empty().set { ch_versions }
        // Cluster the results, varying the resolution ------------------------
        CLUSTER(
            outdir,
            anndata,
            metadata,
            pcs,
            reduced_dims,
            cluster__number_neighbors,
            cluster__methods,
            cluster__resolutions
        )
        ch_versions = ch_versions.mix(CLUSTER.out.versions)

        // Boxplot of phenotype across clusters -------------------------------
        PLOT_PHENOTYPE_ACROSS_CLUSTERS(
            CLUSTER.out.outdir,
            CLUSTER.out.anndata,
            cluster__boxplot_variables
        )
        ch_versions = ch_versions.mix(PLOT_PHENOTYPE_ACROSS_CLUSTERS.out.versions)

        // // Dotplot of marker genes across clusters ----------------------------
        // // Serialize marker genes
        SERIALIZE_KNOWN_MARKERS(
            cluster__known_markers
        )
        PLOT_KNOWN_MARKERS(
            CLUSTER.out.outdir,
            CLUSTER.out.anndata,
            SERIALIZE_KNOWN_MARKERS.out.marker_file
        )
        ch_versions = ch_versions.mix(PLOT_KNOWN_MARKERS.out.versions)
        
        if (params.cluster_validate_resolution_keras){
            CLUSTER_VALIDATE_RESOLUTION_KERAS( 
                CLUSTER.out.outdir,
                CLUSTER.out.anndata,
                CLUSTER.out.metadata,
                CLUSTER.out.pcs,
                CLUSTER.out.reduced_dims,
                CLUSTER.out.clusters,
                cluster_validate_resolution__sparsity,
                cluster_validate_resolution__train_size_cells,
                CLUSTER.out.outdir__reduced_dims
            )
            ch_versions = ch_versions.mix(CLUSTER_VALIDATE_RESOLUTION_KERAS.out.versions)

            PLOT_RESOLUTION_VALIDATE(
                CLUSTER_VALIDATE_RESOLUTION_KERAS.out.plot_input.groupTuple()
            )
            ch_versions = ch_versions.mix(PLOT_RESOLUTION_VALIDATE.out.versions)
        }
            
        // SCCAF(CLUSTER.out.outdir,
        //   CLUSTER.out.anndata,
        //   CLUSTER.out.clusters,
        //   sccaf_minacc)

        // // Generate UMAPs of the results.
        UMAP_CALCULATE_AND_PLOT(
            CLUSTER.out.outdir,
            CLUSTER.out.anndata,
            CLUSTER.out.pcs,
            CLUSTER.out.reduced_dims,
            use_pcs_as_reduced_dims,
            "",
            "cluster",
            n_neighbors,
            umap_init,
            umap_min_dist,
            umap_spread
        )
        ch_versions = ch_versions.mix(UMAP_CALCULATE_AND_PLOT.out.versions)
        dummy_output=UMAP_CALCULATE_AND_PLOT.out.dummy_output
        // // Find marker genes for clusters
        if (params.cluster_markers){
            CLUSTER_MARKERS(
                CLUSTER.out.outdir,
                CLUSTER.out.anndata,
                CLUSTER.out.metadata,
                CLUSTER.out.pcs,
                CLUSTER.out.reduced_dims,
                CLUSTER.out.clusters,
                cluster_marker__methods
            )
            ch_versions = ch_versions.mix(CLUSTER_MARKERS.out.versions)

            // // Find marker genes for clusters using CELLEX
            CELLEX_CLUSTER_MARKERS(
                CLUSTER.out.outdir,
                CLUSTER.out.anndata
            )
            ch_versions = ch_versions.mix(CELLEX_CLUSTER_MARKERS.out.versions)

            // Prep adata file for cellxgene website
            PREP_CELLXGENE(
                CLUSTER.out.outdir,
                CLUSTER.out.anndata
            )
            ch_versions = ch_versions.mix(PREP_CELLXGENE.out.versions)
        }
        emit:
            dummy_output
            versions = ch_versions

}