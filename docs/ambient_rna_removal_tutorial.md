 ## Running and optimising ambient rna removal with cellbender:

You may only want to run the ambient RNA removal and then investigate the reports of the performence of the training.
This is possible with the pipeline by providing these inputs in your params file such as [this one](https://github.com/wtsi-hgi/yascp/blob/v1.5/assets/deploy_scripts/input_setups/cellbender_profile.nf):

```
params {
    do_deconvolution = false
    celltype_assignment.run_celltype_assignment = false
    skip_qc = true
    skip_handover = true
    skip_merge = true
}
```


If after running cellbender you notice that training eppochs were too low, low counts treshold was too high (comon for ATAC and snRNAseq datasets) and you want to tweak some of the samples parameters you can do this by providing the samples params you want to use:
Have a look at the default [params used in the pipeline here](https://github.com/wtsi-hgi/yascp/blob/v1.5/conf/cellbender.conf)
```
params {
    cellbender_rb{
        description = 'Parameters for cellbender remove background.'
        per_sample_thresholds = [
                [ name : 'cellranger700_multi_850906bde9153135a2abd77d0227353e', low_count_threshold : 5, epochs : 100 , learning_rate : 0.0001, zdim : 100, zlayers : 500],
                [ name : 'sample2', low_count_threshold : "", epochs : "" , learning_rate : "", zdim : 100, zlayers : 500],
        ]
}
}
```


<details markdown="1">
<summary><b>Sanger Specific Exacution:</b></summary>

* In Sanger you do not need to set up anything. All you need is an input file:
  ```
      module load HGI/pipelines/yascp/1.5
      yascp cellbender -c input.nf
  ```
