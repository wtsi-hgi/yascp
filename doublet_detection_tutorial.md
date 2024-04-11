## Running Only Doublet Assesments:

You can run doublet assessments exclusively, using tools such as Scrublet, scDblFinder, DoubletDecon, DoubletDetection, and DoubletFinder on your matrix files, if you wish to bypass deconvolution and other processing steps. This functionality is supported by the pipeline.

To facilitate this, you must still specify your input.nf file, which should contain all the sample details from the alignment phase.

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
</details>

By default yascp runs with all of these doublet detection methods, you can switch some of them off by providing these params settings:

```console
    filter_multiplets{
        run_process = true
        doubletDetection{
             run_process = false           
        }
        doubletDecon{
            run_process = false
        }
        scDblFinder{
            run_process = false
        }
        scds{
            run_process = false
        }
        doubletFinder{
            run_process = false
        }
        scrublet{
            run_process= true
        }
    }
```