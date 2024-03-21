 ## Running Only Doublet Assesments:

It is possible to run only the doublet assesments with scrublet, scDblFinder, DoubletDecon, DoubletDetection and DoubletFinder of your matrix files as you may not want to run deconvolution and all other processing steps.
This is possible with the pipeline.
To do this you still have to point to your input.nf file that contains all the samples from allignment

```console
params{
    input_data_table = "/path/to/my/inputs/tsv/input.tsv"
}
```
And then you can run the pipeline as:
```
    nextflow run /path/to/cloned/yascp -profile sanger -entry JUST_DOUBLETS -c input.nf
```

<details markdown="1">
<summary><b>Sanger Specific Exacution:</b></summary>

* In Sanger you do not need to set up anything. All you need is an input file:
  ```
      module load HGI/pipelines/yascp/1.5
      yascp doublets -c input.nf
  ```
