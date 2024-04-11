 ## Cleaning up the intermediate files:

Avoid directly deleting the work directory, as this will result in the loss of all your data. Additionally, be aware that executing this process will clear your cache. Should you need to rerun the pipeline, you will be required to initiate from one of the designated entry points detailed in the tutorials provided, or start anew from the beginning. It is advisable to clean the work directories only after the entire pipeline has been executed and you have confirmed the accuracy of all data, as well as the completion of all processes.

For your convenience, we have integrated an entry point within the pipeline to facilitate the cleaning of the results directory. You can activate this feature by following the instructions provided:
```
nextflow run /path/to/cloned/yascp -profile sanger -c inputs.nf -entry WORK_DIR_REMOVAL
```

The input file is the same file as you used for running the pipeline, no modifications needs to be made.

<details markdown="1">
<summary><b>Sanger Specific Exacution:</b></summary>

* In Sanger you do not need to set up anything. All you need is an input file:
  ```
      module load HGI/pipelines/yascp/1.5
      yascp clean -c input.nf
  ```
</details>