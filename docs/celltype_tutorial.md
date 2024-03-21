 ## Running Only Celltype Assesments:

It is possible to run only the celltype assignemt of your h5ad file as you may not want to run deconvolution and all other processing steps.
This is possible with the pipeline.
To do this you can modify the input.nf file as:

```console
params{
    celltype_assignment{
        run_celltype_assignment=true
        run_keras=true
        run_celltype_assignment=false
        run_azimuth=true
        run_celltypist=false
    }
    skip_preprocessing=true
    file__anndata_merged = '/lustre/scratch126/humgen/projects/sc-eqtl-ibd/analysis/yascp_analysis/2023_08_22-sequencing_pilots/ti_Element/results/merged_h5ad/pre_QC_adata.h5ad'
}
```
And then you can run the pipeline as:
```
    nextflow run /path/to/cloned/yascp -profile sanger -entry JUST_CELLTYPES -c input.nf
```

Aditionally if you have your own celltypist models that you want to use you can edit the default params:
Please take a look on the available models in default [config file](https://github.com/wtsi-hgi/yascp/blob/c55fcfb1a11045e16125f31c20ebe57e0fe81149/conf/qc.conf#L44-L56)
```
params{
    celltypist {
        models = ['/path/to/my/model/Immune_All_High.pkl']
    }
}
```
<details markdown="1">
<summary><b>Sanger Specific Exacution:</b></summary>

* In Sanger you do not need to set up anything. All you need is an input file:
  ```
      module load HGI/pipelines/yascp/1.4
      yascp celltype -c input.nf
  ```
