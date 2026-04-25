// Modules to include.
include {
    PLOT_FILTERED_CELLS;
    PLOT_PREDICTED_SEX;
    PLOT_QC;
    PLOT_DISTRIBUTIONS;PLOT_PCS;
} from "./functions.nf"


workflow PLOT_STATS {
    take:
        file__anndata_merged
        file__cells_filtered
        outdir
        anndata
        n_pcs
        
    main:
        Channel.empty().set { ch_versions }
        // Plot the filtered cells per sample.
        PLOT_FILTERED_CELLS(
            params.outdir,
            file__cells_filtered
        )
        ch_versions = ch_versions.mix(PLOT_FILTERED_CELLS.out.versions)

        // Predict sex from gene expression and check against phenotypes.
        PLOT_PREDICTED_SEX(
            params.outdir,
            file__anndata_merged
        )
        ch_versions = ch_versions.mix(PLOT_PREDICTED_SEX.out.versions)

        // Make QC plots of the merged data.
        PLOT_QC(
            params.outdir,
            file__anndata_merged,
            params.plots_qc.facet_columns.value
        )
        ch_versions = ch_versions.mix(PLOT_QC.out.versions)

        // Plot Distributions
        PLOT_DISTRIBUTIONS(
            params.outdir,
            file__anndata_merged,
            params.plots_qc.facet_columns.value,
            params.plots_qc.variable_columns_distribution_plots.value
        )
        ch_versions = ch_versions.mix(PLOT_DISTRIBUTIONS.out.versions)

        PLOT_PCS(
            params.outdir,
            anndata,
            n_pcs,
            params.umap.colors_quantitative.value,
            params.umap.colors_categorical.value
        )
        ch_versions = ch_versions.mix(PLOT_PCS.out.versions)
        LI = PLOT_PCS.out.out_png
        emit:
            LI
            versions = ch_versions
}

