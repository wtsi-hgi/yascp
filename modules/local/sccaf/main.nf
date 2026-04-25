
include {
    SCCAF_ASSESS_CLUSTERING;
    SCCAF_OPTIMIZE_CLUSTERING;
} from "./functions.nf"

workflow SCCAF {
  take:
    outdir
    anndata
    external_clustering
    min_accuracy

  main:
    Channel.empty().set { ch_versions }
    SCCAF_ASSESS_CLUSTERING(
        outdir,
        anndata,
        external_clustering
    )
    ch_versions = ch_versions.mix(SCCAF_ASSESS_CLUSTERING.out.versions)

    SCCAF_OPTIMIZE_CLUSTERING(
        outdir,
        anndata,
        external_clustering,
        min_accuracy
    )
    ch_versions = ch_versions.mix(SCCAF_OPTIMIZE_CLUSTERING.out.versions)
    emit:
        versions = ch_versions
}
