# nf-core/yascp: Usage
<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->
## Installation
### General Installation Instruction:

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`)

2. Install [`Docker`](https://docs.docker.com/engine/installation/) or [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/) for full pipeline reproducibility.

3. Download/clone the pipeline:

```console
    git clone https://github.com/wtsi-hgi/yascp.git
```
4. The YASCP pipeline is ready to run.

<details markdown="1">
<summary><b>Sanger-Specific Installation:</b></summary>

You don't need to install anything. YASCP is already installed on Farm and can be loaded as a module
  
</details>

## Running the pipeline
### General Execution Instruction:

To run the whole pipeline use the next commands:

For a test dataset run:
```console
   nextflow run /path/to/cloned/yascp -profile test,<docker/singularity,institute>
```
For your dataset run:
```console
   nextflow run /path/to/cloned/yascp -profile <docker/singularity,institute> -c inputs.nf -resume
```

<details markdown="2">
<summary><b>Sanger-Specific Execution:</b></summary>

* If you work on Farm you can run YASCP using the next commands:
  
  Test dataset run:
    ```
      module load HGI/pipelines/yascp/1.6.1
      yascp test
  ```
  Your own dataset run:
  ```
      module load HGI/pipelines/yascp/1.6.1
      yascp -c input.nf
  ```

</details>

## Core Nextflow arguments
To run YASCP you need to specify several core Nextflow arguments like in the example commands above.

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different computing environments.

Several generic profiles are included with the pipeline, guiding it to use software packaged via various methods such as Docker and Singularity (see below). When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g
 <!-- [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/). -->

> You will need to use Docker or Singularity containers for full pipeline reproducibility as currently, we do not support Conda.

<!-- For us it doesn't - but would be nice to do this too ---- The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation). -->

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
    * A generic configuration profile to be used with [Docker](https://docker.com/)
* `singularity`
    * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters
* `institute`
    * A profile with a complete configuration for your institute resources
    * To use your institution profile, replace `institute` with your institution profile name. Many institutions provide profiles (look for yours https://github.com/nf-core/configs/tree/master/conf)
    * If there is no profile for your institution you can create your configure file and specify it using `-c`

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.
You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a config file (including the input declaration config file). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Input declaration config file

This file specifies all inputs to the pipeline and general pipeline parameters. You can find an example input declaration [here](../sample_input/inputs.nf).

Core required/optional inputs are described below. 

```console
params {
    //REQUIRED
    input_data_table = '/path/to/input.tsv' //A samplesheet file containing paths to all the cellranger and pool definition files

    //OPTIONAL
    input = 'cellbender' //This parameter defines whether the ambient RNA removal is skipped ('cellranger') or not ('cellbender'). The default value is 'cellbender'.

    //cellbender_location='/path/to/existing/folder/nf-preprocessing/cellbender' //!!!!! Uncomment this and edit the path, if cellbender results are already available (even partial results). The pipeline will skip the cellbender step for samples that already have cellbender results.

    existing_cellsnp="" //Provide a path to cellsnp results (if they are already available, even partial results) to skip cellsnp step for the files with results.

    genotype_input {
        run_with_genotype_input=true //This parameter defines whether the genotype_input is used (true) or not(false). If this is set to true tsv_donor_panel_vcfs has to be specified
        tsv_donor_panel_vcfs = "/path/to/reference/panel/vcf_inputs.tsv" //A file containing paths to vcf files with a priori known genotypes that we want to compare the genotypes from samples with
    }
}


```
### Required parameters
`input_data_table` - a samplesheet file containing paths to all the cellranger and pool definition files.

### Optional parameters
`input` - This parameter defines whether the ambient RNA removal is skipped ('cellranger') or not ('cellbender'). The default value is 'cellbender'. This option can be useful if you can't use GPUs. For more details see `Tips to avoid rerunning time-consuming parts of the pipeline.` below.

`cellbender_location` - uncomment this and edit the path, if cellbender results are already available (even partial results). The pipeline will skip the cellbender step for samples that already have cellbender results. For more details see `Tips to avoid rerunning time-consuming parts of the pipeline.` below.

`existing_cellsnp` - provide a path to cellsnp results (if they are already available, even partial results) to skip cellsnp step for the files whth results. For more details see `Tips to avoid rerunning time-consuming parts of the pipeline.` below.

`run_with_genotype_input` - this parameter defines whether the genotype_input is used (true) or not(false). If this is set to true `tsv_donor_panel_vcfs` has to be specified.

`tsv_donor_panel_vcfs` - a file containing paths to vcf files with a priori known genotypes that we want to compare the genotypes from samples with.

### You can find more optional parameters [here](optional_parameters.md).

## Samplesheet input
This file specifies sample IDs, the number of pooled donors, IDs of individuals with priori known genotypes, and paths to 10x files.
It has to be a tab-separated file with 4 columns and a header as shown in the example below.
You can find an example samplesheet [here](../sample_input/input_table.tsv).


| experiment_id   | n_pooled | donor_vcf_ids    |  data_path_10x_format   |
|-----------------|----------|------------------|-------------------------|
| Pool1 |   1      | "id3"            | path/to/10x_folder      |
| Pool2|   2      | "id1,id2"        | path/to/10x_folder      |
 

Columns description:
* **experiment_id** - is the name of the sample
* **n_pooled** - indicates how many donors are pooled in the 10x run (if only 1 then scrubblet will be used to remove doublets)
* **donor_vcf_ids** - if using genotypes, here an id of individuals can be added to subset VCFs used to deconvolute samples (need to be as listed in VCF file provided)
* **data_path_10x_format** - path to a 10x folder containing bam, bai, metrics_summary.csv files and raw_barcodes folder

**path/to/10x_folder** can contain output files from both cellranger 6 and cellranger 7. Overall we need the following files for the pipeline to run smoothly:

```console
10x_folder/
    ./possorted_genome_bam.bai
    ./possorted_genome_bam.bam
    ./raw_feature_bc_matrix
        ./matrix.mtx.gz
        ./features.tsv.gz
        ./barcodes.tsv.gz
    ./filtered_feature_bc_matrix
        ./matrix.mtx.gz
        ./features.tsv.gz
        ./barcodes.tsv.gz
    ./metrics_summary.csv
    ./web_summary.html
    ./molecule_info.h5
```
You could also provide a path to this file by using a flag:
```console
--input_data_table '[path to samplesheet file]'
```

## Genotypesheet input (optional)
This file contains paths to VCFs and cohort labels associated with them.
A genotypesheet can be provided to the pipeline to improve sample deconvolution and detect whether the sample you have is the sample you are expecting (through GT matching).
The pipeline will determine which cohort the deconvoluted sample comes from (if any).

In the following example, we have 3 cohorts: Cohort1 has genotypes for each of the chromosomes - this is acceptable, as the pipeline will use all chromosome files to identify whether the sample is part of this cohort. The other 2 cohorts have a merged VCF file for all the chromosomes. This is also acceptable, as it will determine whether the sample belongs to this cohort in one step. After evaluating all cohorts the pipeline will assign the sample to the single donor that is the most likely real match.

You can find an example genotypesheet [here](../sample_input/vcf_inputs.tsv).

| label   | vcf_file_path    |
|-----------------|----------|
| Cohort1 |   /ful/path/to/vcf_bcf/file/in/hg38/format/without/chr/prefix/chr1.vcf.gz      |
| Cohort1 |   /ful/path/to/vcf_bcf/file/in/hg38/format/without/chr/prefix/chr2.vcf.gz      |
| .... |   ....      |
| Cohort2 |   /ful/path/to/vcf_bcf/file/in/hg38/format/without/chr/prefix/full_cohort2_for_all_chr.vcf.gz      |
| Cohort3 |   /ful/path/to/vcf_bcf/file/in/hg38/format/without/chr/prefix/full_cohort2_for_all_chr.vcf.gz      |


## Tips to avoid rerunning time-consuming parts of the pipeline.
To avoid rerunning time-consuming steps of the pipeline you can specify the next parameters in the input declaration config file:

### 1. input
You can skip the cellbender step by adding `input = 'cellranger'` to the input declaration config file. You might consider this option because the cellbender step is time-consuming and requires GPUs.
The pipeline will skip ambient RNA removal and proceed with deconvolution based on cellranger. For more details see [optional parameters](Optional_parameters.md)

``` console
params{
    input_data_table = '/path/to/input.tsv' //A samplesheet file containing paths to all the cellranger and pool definition files
    input = 'cellranger'
}
```

### 2. cellbender_location
You can avoid running cellbender multiple times if you have complete or partial cellbender results.
If you specify a path to the folder with cellbender results in the input declaration config file, cellbender will be run on all the samples without results.
```
params{
    input_data_table = '/path/to/input.tsv' //A samplesheet file containing paths to all the cellranger and pool definition files
    cellbender_location='/full/path/to/results/nf-preprocessing/cellbender'
}
```
The cellbender results folder structure should look like this: 
```console
    Sample1
    Sample2
    Sample3
    qc_cluster_input_files
        file_paths_10x-*FPR_0pt1
        file_paths_10x-*FPR_0pt05
        file_paths_10x-*FPR_0pt01
```

### 3. existing_cellsnp
You can avoid running cellsnp multiple times if you have complete or partial cellsnp results.
If you specify a path to cellsnp files in the input declaration config file, cellsnp will work only with the files that haven't yet been run:

``` console
params{
    input_data_table = '/path/to/input.tsv' //A samplesheet file containing paths to all the cellranger and pool definition files
    existing_cellsnp='/full/path/to/results/cellsnp'
}
```
<!-- 2. full_vcf_file = points to vcf file to be used.
4. subset_genotypes = indicates to subset genotypes for an input to be used in Vireo.
5. run_celltype_assignment = runs celltypist and Azimuth if PBMC data is used.
6. file__anndata_merged = if the preprocessing is already done it can input a merged h5ad which will proceed directly to integration and clustering.
7. extra_metadata = any extra metadata to be added for samples.
8. input_data_table = is a file pointing to the 10x files as per: -->

## Pipeline custom configuration
If you need to customise the pipeline please read **[Custom configuration](custom_configuration.md)** for more details.

## Reproducibility
It is a good idea to specify a pipeline version (or a checkout tag indicated when running `git log`) when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.
<!-- TODO - add a description about reproducibility something like this: currently we don't have a release;
It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/yascp releases page](https://github.com/nf-core/yascp/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline so that you'll know what you used when you look back in the future. -->

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or a similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted by your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```
