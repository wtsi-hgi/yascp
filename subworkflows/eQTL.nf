
include {PSEUDOBULK_AGREGATION} from "$projectDir/modules/nf-core/modules/pseudobulk_aggregation/main"

workflow eQTL{
    take:
        file__anndata_merged
        assignments_all_pools
    main:
        PSEUDOBULK_AGREGATION(file__anndata_merged,params.eQTL.aggregation_collumn,params.eQTL.n_min_cells,params.eQTL.n_min_individ,assignments_all_pools)
        // Here we need to pass in the GT match results if the Genotypes are used, otherwise dont do anything.

}