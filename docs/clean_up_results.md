 ## Cleaning up the intermediate files:

Please dont just delete the work directory - all your data will be lost. 
Also note that after running this process your cache will be lost and you in case you need to rerun pipeline you will either have to start from one of the entry points as described in these tutorials or start from scratch.
After a full pipeline is run and you have validated that all the data is accurate and alll the processes finished then its a good idea to clean up the work directories.

We have provided an entry point in pipeline to clean up the results directory. 
You can run this by:
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
