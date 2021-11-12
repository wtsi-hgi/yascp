
include {
    sccaf_assess_clustering;
    sccaf_optimize_clustering;
} from "./functions.nf"

workflow SCCAF {
  take:
    outdir
    anndata
    external_clustering
    min_accuracy

  main:
   
    sccaf_assess_clustering(
        outdir,
        anndata,
        external_clustering
    )

    sccaf_optimize_clustering(
        outdir,
        anndata,
        external_clustering,
        min_accuracy
    )
}
