## Running Only Celltype Assignments or Adjusting params for a pipeline regarding celltype assignments:

You have the flexibility to execute solely the cell type assignment for your h5ad file, bypassing the need for deconvolution and other processing steps. This customization is achievable within the pipeline framework. To tailor the pipeline to perform only the cell type assignment, you can adjust the input.nf file as follows:

```console
params{
    celltype_assignment{
        run_celltype_assignment=true
        run_keras=true
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

<details markdown="1">
<summary><b>Sanger Specific Exacution:</b></summary>

* In Sanger you do not need to set up anything. All you need is an input file:
  ```
      module load HGI/pipelines/yascp/1.7
      yascp celltype -c input.nf
  ```
</details>

### Celltypist
Additionally, if you possess custom CellTypist models you wish to utilize, you can modify the default parameters. For an overview of the models currently available by default, please refer to the [config file](https://github.com/wtsi-hgi/yascp/blob/c55fcfb1a11045e16125f31c20ebe57e0fe81149/conf/qc.conf#L44-L56)
```
params{
    celltypist {
        models = ['/path/to/my/model/Immune_All_High.pkl']
    }
}
```

### Azimuth 
Additionally, if you have a custom Azimuth model or one obtained from [Azimuth Zenodo](https://azimuth.hubmapconsortium.org/references/)  that you wish to use, you can modify the default parameters. By default, the pipeline uses a PBMC reference, but you have the flexibility to change this and include any references that meet your requirements.
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

