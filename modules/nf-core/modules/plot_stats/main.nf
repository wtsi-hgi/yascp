// Modules to include.
include {
    plot_filtered_cells;
    plot_predicted_sex;
    plot_qc;
    plot_distributions;plot_pcs;
} from "./functions.nf"


workflow PLOT_STATS {
    take:
        file__anndata_merged
        file__cells_filtered
        outdir
        anndata
        n_pcs
        
    main:
        // Plot the filtered cells per sample.
        plot_filtered_cells(
            params.outdir,
            file__cells_filtered
        )

        // Predict sex from gene expression and check against phenotypes.
        plot_predicted_sex(
            params.outdir,
            file__anndata_merged
        )

        // Make QC plots of the merged data.
        plot_qc(
            params.outdir,
            file__anndata_merged,
            params.plots_qc.facet_columns.value
        )

        // Plot Distributions
        plot_distributions(
            params.outdir,
            file__anndata_merged,
            params.plots_qc.facet_columns.value,
            params.plots_qc.variable_columns_distribution_plots.value
        )

        plot_pcs(
            outdir,
            anndata,
            n_pcs,
            params.umap.colors_quantitative.value,
            params.umap.colors_categorical.value
        )
        LI = plot_pcs.out.out_png
        emit:
            LI
}

