# nf-core/yascp: Usage
## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->
## Installation
YASC requires 
[Nextflow](https://www.nextflow.io/)

and one of the container platforms
[Docker](https://docker.com/) or [Singularity](https://sylabs.io/docs/)

You only need to clone the repository:
```console
    git clone https://github.com/wtsi-hgi/yascp.git
```
The YASCP pipeline is ready to run.

## Running the pipeline

<details markdown="1">
<summary><b>Sanger-Specific Execution (you don't need to set up anything):</b></summary>


  Test dataset run:
  ```
      module load HGI/pipelines/yascp/1.6.1
      yascp test
  ```
  Your dataset run:
  ```
      module load HGI/pipelines/yascp/1.6.1
      yascp -c input.nf
  ```

</details>

To run the whole pipeline use the next commands:

For a test dataset run:
```console
   nextflow run /path/to/cloned/yascp -profile sanger,test,singularity
```
For your dataset run:
```console
   nextflow run /path/to/cloned/yascp -profile sanger,singularity -c inputs.nf -resume
```

These commands will launch the pipeline with the `sanger` configuration profile. Instead of the sanger profile please select your institution profile - the available profiles are here: https://github.com/nf-core/configs/tree/master/conf .
More information about profiles is bellow.

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different computing environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g
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

```#COMMENT what about institute profiles?```

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.
You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Input declaration config file
An [example samplesheet](../sample_input/inputs.nf) has been provided with the pipeline.

Since we have multiple inputs in the pipeline we point to each of them in a sample config file. There are multiple required/optional inputs that are described below.

```console
params {
    extra_metadata = '/path/to/extra_metadata.tsv'   //Sometimes users may want to merge extra known metadata for a pool in the h5ad files prior to qc

    extra_sample_metadata ="/path/to/donor_extra_metadata.tsv"  //Sometimes users may want to merge extra known metadata for a donor within pool prior to qc

    input_data_table = '/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/OneK1k/onek1k_test_dataset/input.tsv' //Required!! This points to all the cellranger files and pool definition files.

    split_ad_per_bach=true //if not splitting the celltype assignment will be run on full dataset together
    //cellbender_location='/path/to/existing/folder/nf-preprocessing/cellbender' //!!!!! uncoment and change path if already have results - if cellbender is run already then can skip this by selecting  input = 'existing_cellbender' instead input = 'cellbender'

    existing_cellsnp="" // if we have run cellsnp before we can skip this process by letting yascp capture the files

    genotype_input {
        run_with_genotype_input=true //if false do not need the genotype_input parameters.
        vireo_with_gt=false // Vireo is capable in runing both with genotypes and without. Here we define in which mode we want to run it.
        posterior_assignment = false //if this is set to true, we will perform the genotype donor matching after the deconvolution is performed.
        subset_genotypes = false
        tsv_donor_panel_vcfs = "/path/to/reference/panel/vcf_inputs.tsv" //this is a panel of vcf files that we want to compar the genotypes with
    }
}


```

This input.nf file will be provided when the pipeline is executed:
    ```
    nextflow run /path/to/cloned/nfCore_scRNA -profile sanger -resume -c input.nf
    ```

## Samplesheet input
An [example samplesheet](../sample_input/input_table.tsv) has been provided with the pipeline.
As per above main file required is a file containing paths to 10x files in a format:


| experiment_id   | n_pooled | donor_vcf_ids    |  data_path_10x_format   |
|-----------------|----------|------------------|-------------------------|
| Pool1 |   1      | "id3"            | path/to/10x_folder      |
| Pool2|   2      | "id1,id2"        | path/to/10x_folder      |

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a tab-separated file with 4 columns and a header as shown in the example above.

Where:
* **experiment_id** - is the name of the sample
* **n_pooled** - indicates how many donors are pooled in the 10x run (if only 1 then scrubblet will be used to remove doublets)
* **donor_vcf_ids** - if using genotypes, here an id of individuals can be added to subset VCFs used to deconvolute samples (need to be as listed in VCF file provided)
* **data_path_10x_format** - path to a 10x folder containing bam, bai, metrics_summary.csv files and raw_barcodes folder

**path/to/10x_folder** these outputs can be outputs from both cellranger 6 and cellranger 7. Overall we need the following files for the pipeline to run smoothly:

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
An [example genotypesheet](../sample_input/vcf_inputs.tsv) has been provided with the pipeline.
Genotypesheet can be provided to the pipeline to perform a better sample deconvolution and detect whether the sample you are expecting is really the sample (through the GT match).```#COMMENT I didn't understand the part with "is really the sample"```
The pipeline will figure out which cohort the deconvoluted sample comes from (if any). In the following example, we have 3 cohorts: Cohort1 has genotypes for each of the chromosomes - this is ok as the pipeline will use all chromosome files to figure out whether the sample is part of this. The other 2 cohorts have a merged VCF file for all the chromosomes. This is also ok as it will figure out whether the sample belongs to this cohort in one go. After looking at all these cohorts the pipeline will assign only 1 donor corresponding to which one is the most likely real match

| label   | vcf_file_path    |
|-----------------|----------|
| Cohort1 |   /ful/path/to/vcf_bcf/file/in/hg38/format/without/chr/prefix/chr1.vcf.gz      |
| Cohort1 |   /ful/path/to/vcf_bcf/file/in/hg38/format/without/chr/prefix/chr2.vcf.gz      |
| .... |   ....      |
| Cohort2 |   /ful/path/to/vcf_bcf/file/in/hg38/format/without/chr/prefix/full_cohort2_for_all_chr.vcf.gz      |
| Cohort3 |   /ful/path/to/vcf_bcf/file/in/hg38/format/without/chr/prefix/full_cohort2_for_all_chr.vcf.gz      |


## Extra pool metadata sheet (optional)
An [example pool metadata](../sample_input/extra_metadata.tsv) has been provided with the pipeline.

Users may want to provide extra metadata for each of the pools that can be used for clustering, regression or plotting purposes.

| experiment_id   | Experimental design | Library prep date | Stimulation time    | ...   |
|-----------------|----------|------------------|-------------------------|-----|
| Pool1 |   1      |   20/01/2023          | 24h      |  |
| Pool2|   2      | 21/01/2023        | 48h      |  |

## Extra donor within pool metadata sheet (optional)
An [example metadata for donors in a pool](../sample_input/extra_metadata_donors.tsv) has been provided with the pipeline.

If users have used genotypes in the pipeline then upon deconvolution and gt match we will be able to tell which donor is which. In this case, if users have any extra information for each of the donors within a pool then this extra metadata information can also be provided in the same format as above. To make sure that the correct metadata gets attached to the correct donor the experiment_id should contain experiment_id__donor_genotype_id  
(Note: if you provided a bridging file this should be experiment_id__phenotype_id) 

| experiment_id   | Sex | Age | Condition    | ...   |
|-----------------|----------|------------------|-------------------------|-----|
| Pool1__donor1 |   M      |   67          | PAH      |  |
| Pool1__donor2|   M      | 22        | CD      |  |
| Pool1__donor3|   F      | 43        | CD      |  |
| ...|   ...      | ...        | ...      |  |
| Pool2__donor1|   F      | 12        | PAH      |  |
| ...|   ...      | ...        | ...      | ... |
| Pool2__donorN|   M      | 88        | AH      |  |

## Genotype to phenotype bridging file (optional)
An [genotype to phenotype bridging file](../sample_input/genotype_phenotype_bridge.tsv) has been provided with the pipeline.

Sometimes IDs that we expect in our [input files](../sample_input/genotype_phenotype_bridge.tsv) donor_vcf_ids may correspond to phenotype IDs instead of genotype IDs. Since the pipeline performs the checks of whether the donor that we get is the one we expect according to this field (very important step for the Cardinal project) we want to map the genotype IDs to phenotype IDs. This will be handled by the pipeline.

| oragene_id   | s00046_id    |
|-----------------|----------|
| 682_683 |   pheno_682_683      |
| 684_684 |   pheno_682_683      |
| .... |   ....      |


## Some tricks to avoid rerunning the pipeline over and over if you already have some partial data

1. You can avoid running cellbender multiple times. If you have even partial cellbender results you can provide a path to the folder that contains them. cellbender will be run on all the samples besides the ones that are captured by [cellbender_location='/full/path/to/results/nf-preprocessing/cellbender']. Other options - [cellranger] - which avoids ambient RNA removal and proceeds with deconvolution based on cellranger. If you are providing a path to cellbender_location ='??' - specify location to the results directory containing:
```
params{
    cellbender_location='/full/path/to/results/nf-preprocessing/cellbender'
}
```
This should contain: 
```console
    Sample1
    Sample2
    Sample3
    qc_cluster_input_files
        file_paths_10x-*FPR_0pt1
        file_paths_10x-*FPR_0pt05
        file_paths_10x-*FPR_0pt01
```

2. existing_cellsnp = '' - If you point to the path of partial cellsnp files these will be captured in the pipeline and utilised in downstream processes, and only cellsnp of the files that don't have the runs performed on cellsnp will proceed:

``` console
params{
    existing_cellsnp='/full/path/to/results/cellsnp'
}
```
<!-- 2. full_vcf_file = points to vcf file to be used.
4. subset_genotypes = indicates to subset genotypes for an input to be used in Vireo.
5. run_celltype_assignment = runs celltypist and Azimuth if PBMC data is used.
6. file__anndata_merged = if the preprocessing is already done it can input a merged h5ad which will proceed directly to integration and clustering.
7. extra_metadata = any extra metadata to be added for samples.
8. input_data_table = is a file pointing to the 10x files as per: -->

### Reproducibility
It is a good idea to specify a pipeline version (or a checkout tag indicated when running `git log`) when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.
<!-- TODO - add a description about reproducibility something like this: currently we don't have a release;
It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/yascp releases page](https://github.com/nf-core/yascp/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline so that you'll know what you used when you look back in the future. -->


