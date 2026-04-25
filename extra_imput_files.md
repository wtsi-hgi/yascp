## Extra pool metadata sheet (optional)
This file contains extra metadata for each of the pools that can be used for clustering, regression or plotting purposes.

You can find an example file with pool metadata [here](../sample_input/extra_metadata.tsv).

| experiment_id   | Experimental design | Library prep date | Stimulation time    | ...   |
|-----------------|----------|------------------|-------------------------|-----|
| Pool1 |   1      |   20/01/2023          | 24h      |  |
| Pool2|   2      | 21/01/2023        | 48h      |  |

## Extra donor within pool metadata sheet (optional)
This file contains extra metadata for each of the donors within a pool. To make sure that the correct metadata gets attached to the correct donor the experiment_id should contain experiment_id__donor_genotype_id  

You can find an example file with metadata for donors in a pool [here](../sample_input/extra_metadata_donors.tsv).

(Note: if you provided a bridging file this should be experiment_id__phenotype_id) 

| experiment_id   | Sex | Age | Condition    | ...   |
|-----------------|----------|------------------|-------------------------|-----|
| Pool1__donor1 |   M      |   67          | PAH      |  |
| Pool1__donor2|   M      | 22        | CD      |  |
| Pool1__donor3|   F      | 43        | CD      |  |
| ...|   ...      | ...        | ...      |  |
| Pool2__donor1|   F      | 12        | PAH      |  |
| ...|   ...      | ...        | ...      | ... |
| Pool2__donorN|   M      | 88        | AH      |  |

## Genotype to phenotype bridging file (optional)

Sometimes IDs that we expect in donor_vcf_ids column of our samplesheet may correspond to phenotype IDs instead of genotype IDs. Since the pipeline performs the checks of whether the donor that we get is the one we expect according to this field (very important step for the Cardinal project) we want to map the genotype IDs to phenotype IDs. This will be handled by the pipeline.

You can find an example genotype to phenotype bridging file [here](../sample_input/genotype_phenotype_bridge.tsv).

| oragene_id   | s00046_id    |
|-----------------|----------|
| 682_683 |   pheno_682_683      |
| 684_684 |   pheno_682_683      |
| .... |   ....      |

