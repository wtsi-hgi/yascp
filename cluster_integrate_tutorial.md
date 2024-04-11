## Performing only Integration, Clustering and Cluster Assesments
You may wish to focus solely on integration using BBKNN, Harmony, or Seurat, as well as on clustering and cluster assessments. Alternatively, you might prefer to integrate based on a different variable and/or produce clusters at varying resolutions.

With Yascp, this level of customization is achievable by adjusting the necessary parameters in the settings. There are numerous parameters available; however, you do not need to modify all of them—simply change those that are relevant to your specific requirements:

```console
params{
    skip_preprocessing=true
    file__anndata_merged = '/path/to/merged/andata/pre_QC_adata.h5ad'

    cluster{
        description = """Parameters for clustering. All pairwise combinations of
        method and resolution will be performed."""
        number_neighbors{
            description = """Number of neighbors. If <= 0, uses number of unique
            experiment_id."""
            value = 15
        }
        methods{
            description = 'Clustering method. Valid options [leiden|louvain].'
            value = 'leiden'
        }
        resolutions{
            description = 'Clustering resolution.'
            value = [0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
        }

        variables_boxplot{
            decription = 'Generate boxplots of these variables for each cluster.'
            value ='n_cells,total_counts,pct_counts_gene_group__mito_transcript'
        }
    }


    umap{
        run_process = true
        colors_quantitative{
            description = 'Comma separated string of quantitative variables that will be used to color points.'
            value = 'n_cells,total_counts,pct_counts_gene_group__mito_transcript,prob_doublet,pct_counts_gene_group__ribo_rna,Azimuth:predicted.celltype.l2.score,Azimuth:mapping.score,log10_ngenes_by_count'
        }
        colors_categorical{
            description = 'Comma separated string of categorical variables that will be used to color points.'
            value = 'cell_passes_qc,cell_passes_qc-per:Azimuth:L0_predicted.celltype.l2,experiment_id,Azimuth:predicted.celltype.l2,Celltypist:Immune_All_Low:predicted_labels,Celltypist:Immune_All_High:predicted_labels,donor_id'
        }
    }

    harmony{
        run_process= true
        description= 'Parameters for harmony.'
        variables_and_thetas{
            description = 'Tuples of metadata columns and corresponding thetas.'
            value = [
                [ variable: 'experiment_id', theta: 1.0 ],
                ]
        }
    }

    bbknn{
        run_process = true
        description = 'Parameters for BBKNN'
        batch_variable{
            description= 'Variable to use for batch correction'
            value = 'experiment_id'
        }
    }

}
```
If you have already completed normalization and outlier filtering and simply wish to proceed with clustering, you can specify this by adjusting the corresponding parameter. Here's how you can provide this setting:
```
params{
    dont_integrate_just_cluster = true
}
```

<details markdown="1">
<summary><b>Sanger Specific Exacution:</b></summary>

* In Sanger you do not need to set up anything. All you need is an input file:
  ```
      module load HGI/pipelines/yascp/1.5
      yascp clustering -c input.nf
  ```

