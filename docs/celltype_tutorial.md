 ## Running Only Celltype Assignments or Adjusting params for a pipeline regarding celltype assignments:

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

### Celltypist
Aditionally if you have your own celltypist models that you want to use you can edit the default params:
Please take a look on the available models in default [config file](https://github.com/wtsi-hgi/yascp/blob/c55fcfb1a11045e16125f31c20ebe57e0fe81149/conf/qc.conf#L44-L56)
```
params{
    celltypist {
        models = ['/path/to/my/model/Immune_All_High.pkl']
    }
}
```

### Azimuth 
Aditionally if you have your own azimuth model or a model retrieved from [Azimuth Zenodo](https://azimuth.hubmapconsortium.org/references/) that you want to use you can edit the default params. By default we run a PBMC reference in pipeline but you can change this and add any references you are interested in.
```
params{
    azimuth{
        run_process = true
        celltype_refsets = [
                //# [ name : 'kidney', refset : "/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/jaguar_yascp/nieks_pipeline/yascp_run/ref_kidney", annotation_labels : "cluster,subclass" ],
                [ name : 'PBMC', refset : "PBMC", annotation_labels : "celltype.l2,celltype.l1,celltype.l3" ],
            ]
    }
}
```

<details markdown="1">
<summary><b>Sanger Specific Exacution:</b></summary>

* In Sanger you do not need to set up anything. All you need is an input file:
  ```
      module load HGI/pipelines/yascp/1.5
      yascp celltype -c input.nf
  ```
