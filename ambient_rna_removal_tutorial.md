## Running and optimising ambient rna removal with cellbender:

If your goal is to execute just the ambient RNA removal step and subsequently analyze the performance reports of the training, the pipeline accommodates this process. By specifying the appropriate inputs in your parameters file, similar to the example found in [this one](https://github.com/wtsi-hgi/yascp/blob/v1.5/assets/deploy_scripts/input_setups/cellbender_profile.nf) you can tailor the pipeline to focus solely on ambient RNA removal and performance evaluation.

```
params {
    do_deconvolution = false
    celltype_assignment.run_celltype_assignment = false
    skip_qc = true
    skip_handover = true
    skip_merge = true
}
```


If you observe that the training epochs were insufficiently low, or that the count threshold was too high (a common issue with ATAC and snRNAseq datasets) after running CellBender, you can adjust the sample parameters as needed. To customize the parameters for your samples, you can refer to the default settings used in the pipeline, which are available in the [params used in the pipeline here](https://github.com/wtsi-hgi/yascp/blob/v1.5/conf/cellbender.conf)
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

So your final input.nf file will look something like this:
```
params {
    input_data_table = '/path/to/input.tsv' //Required!! This points to all the cellranger files and pool definition files.
    do_deconvolution = false
    celltype_assignment.run_celltype_assignment = false
    skip_qc = true
    skip_handover = true
    skip_merge = true

    cellbender_rb{
        description = 'Parameters for cellbender remove background.'
        per_sample_thresholds = [
                [ name : 'cellranger700_multi_850906bde9153135a2abd77d0227353e', low_count_threshold : 5, epochs : 100 , learning_rate : 0.0001, zdim : 100, zlayers : 500],
                [ name : 'sample2', low_count_threshold : "", epochs : "" , learning_rate : "", zdim : 100, zlayers : 500],
        ]
    }
}
```
And then you can run the pipeline as:
```
    nextflow run /path/to/cloned/yascp -profile sanger -entry JUST_CELLBENDER -c input.nf
```

<details markdown="1">
<summary><b>Sanger Specific Exacution:</b></summary>

* In Sanger you do not need to clone the repo. All you need is an input file:
  ```
      module load HGI/pipelines/yascp/1.5
      yascp cellbender -c input.nf
  ```
