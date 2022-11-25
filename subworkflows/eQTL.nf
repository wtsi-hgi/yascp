
include {PSEUDOBULK_AGREGATION} from "$projectDir/modules/nf-core/modules/pseudobulk_aggregation/main"

workflow eQTL{
    take:
        file__anndata_merged
    main:
        PSEUDOBULK_AGREGATION(file__anndata_merged,params.eQTL.aggregation_collumn,params.eQTL.n_min_cells,params.eQTL.n_min_individ)
}