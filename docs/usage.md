# nf-core/yascp: Usage
## Introduction

<!-- TODO nf-core: Add documentation about anything specific to running your pipeline. For general topics, please point to (and add to) the main nf-core website. -->

## Input declaration config file
An [example samplesheet](../sample_input/inputs.nf) has been provided with the pipeline.

Since we have multiple inputs in pipeline we point to each of them in a sample config file. There are multiple required/optional inputs that are described bellow.

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

This file will be provided when pipeline is executed:
    ```
    nextflow run /path/to/cloned/nfCore_scRNA -profile sanger -resume -c input.nf
    ```

## Samplesheet input
An [example samplesheet](../sample_input/input_table.tsv) has been provided with the pipeline.
As per above main file required is a paths to 10x files in a format:


| experiment_id   | n_pooled | donor_vcf_ids    |  data_path_10x_format   |
|-----------------|----------|------------------|-------------------------|
| Pool1 |   1      | "id3"            | path/to/10x_folder      |
| Pool2|   2      | "id1,id2"        | path/to/10x_folder      |

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use this parameter to specify its location. It has to be a tab-separated file with 3 columns, and a header row as shown in the examples below.

Where:
* **experiment_id** - is the name of the sample
* **n_pooled** - indicates how many donors are pooled in the 10x run (if only 1 then scrubblet will be used to remove doublets)
* **donor_vcf_ids** - if using genotyes, here an id of individuals can be added to subset vcfs used to deconvolute samples (need to be as listed in vcf file provided)
data_path_10x_format - path to a 10x folder containing bam, bai, metrics_summary.csv files and raw_barcodes folder

**path/to/10x_folder** these outputs can be an output froom both cellranger 6 and cellranger 7. Overall we need the folowing files for pipeline to run smoothly:

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
You could also provide path to this file by providing a flag:
```console
--input_data_table '[path to samplesheet file]'
```

## Genotypesheet input (optional)
An [example genotypesheet](../sample_input/vcf_inputs.tsv) has been provided with the pipeline.
Genotypesheet can be provided to pipeline to perform a beter sample deconvolution, detect wheather the sample you are expecting is really the sample (through the GT match).
Pipeline will figure out which cohort the deconvoluted sample comes from (if any). In the folowing exapme we have 3 cohorts: Cohort1 has genotypes for each of the chromosomes - this is ok as pipeline will use all chromosome files to figure out whether the sample is part of this. Other 2 cohorts has a merged vcf file for all the chromosomes. This is alos ok as it will figure out whether the sample bellongs to this cohort in one go. After looking at all these cohort pipeline will asign only 1 donor corresponding to which one is the most likely real match

| label   | vcf_file_path    |
|-----------------|----------|
| Cohort1 |   /ful/path/to/vcf_bcf/file/in/hg38/format/without/chr/prefix/chr1.vcf.gz      |
| Cohort1 |   /ful/path/to/vcf_bcf/file/in/hg38/format/without/chr/prefix/chr2.vcf.gz      |
| .... |   ....      |
| Cohort2 |   /ful/path/to/vcf_bcf/file/in/hg38/format/without/chr/prefix/full_cohort2_for_all_chr.vcf.gz      |
| Cohort3 |   /ful/path/to/vcf_bcf/file/in/hg38/format/without/chr/prefix/full_cohort2_for_all_chr.vcf.gz      |


## Extra pool metadata sheet (optional)
An [example pool metadata](../sample_input/extra_metadata.tsv) has been provided with the pipeline.

Users may want to provide extra metadata for each of the pools that can be used in clustering, regression or for plotting purposes.

| experiment_id   | Experimental design | Library prep date | Stimulation time    | ...   |
|-----------------|----------|------------------|-------------------------|-----|
| Pool1 |   1      |   20/01/2023          | 24h      |  |
| Pool2|   2      | 21/01/2023        | 48h      |  |

## Extra donor within pool metadata sheet (optional)
An [example metadata for donors in pool](../sample_input/extra_metadata_donors.tsv) has been provided with the pipeline.

If users have used genotypes in pipeline then upon deconvolution and gt match we will be able to tell which donor is which. In this case if users have any extra information for each of the donors within pool then this extra metadata information can also be provided in same format as above. To make sure that the correct metadata gets attached to the correct donor the experiment_id should contain experiment_id__donor_genotype_id  
(Note: if you provided bridging file this should be )experiment_id__phenotype_id 

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

Sometimes IDs that we expect in our [input files](../sample_input/genotype_phenotype_bridge.tsv) donor_vcf_ids may corespond to phenotype IDs instead of genotype ids. Since pipeline performs the checks of whether the donor that we get is the one we expect acording to this field (very important step for Cardinal project) we want to map the genotype ids to phenotype ids. This will be handled by the pipeline.

| oragene_id   | s00046_id    |
|-----------------|----------|
| 682_683 |   pheno_682_683      |
| 684_684 |   pheno_682_683      |
| .... |   ....      |


## Some tricks to avoid reruning pipeline over and over if you already have some partial data

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

2. existing_cellsnp = '' - If you point to the path of partial cellsnp files these will be captured in pipeline and utilised in downstram processes, and only cellsnp of the files that dont have the runs performed on cellsnp will proceed:

params{
    existing_cellsnp='/full/path/to/results/cellsnp'
}

<!-- 2. full_vcf_file = points to vcf file to be used.
4. subset_genotypes = indicates to subset genotypes for an input to be used in Vireo.
5. run_celltype_assignment = runs celltypist and Azimuth if PBMC data is used.
6. file__anndata_merged = if all preprocession has already been doe can input a marged h5ad which will skio all the cellbender and deconvolution.
7. extra_metadata = any extra metadata to be added for samples.
8. input_data_table = is a file pointing to the 10x files as per: -->


## Running the pipeline

The typical command for running the pipeline is as follows [Instead of sanger please select your institution] - the available profiles are here: https://github.com/nf-core/configs/tree/master/conf . you can also provide your own config file by adding an extra input file -c /path/to/my/default/config/file.nf:

```console
nextflow run yascp -profile sanger -c inputs.nf --nf_ci_loc $PWD -resume > nextflow.nohup.log 2>&1 & 
```

This will launch the pipeline with the `sanger` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```console
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

### Reproducibility
It is a good idea to specify a pipeline version (or a checkout tag indicated when runing `git log`) when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.
<!-- TODO - add a description about reproducability something like this: currently we dont have a release;
It is a good idea to specify a pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [nf-core/yascp releases page](https://github.com/nf-core/yascp/releases) and find the latest version number - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future. -->

## Core Nextflow arguments

> **NB:** These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen).

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Conda) - see below. When using Biocontainers, most of these software packaging methods pull Docker containers from quay.io e.g
 <!-- [FastQC](https://quay.io/repository/biocontainers/fastqc) except for Singularity which directly downloads Singularity images via https hosted by the [Galaxy project](https://depot.galaxyproject.org/singularity/) and Conda which downloads and installs software locally from [Bioconda](https://bioconda.github.io/). -->

> You will need to use Docker or Singularity containers for full pipeline reproducibilityas currently we do not support Conda.

<!-- For us it doesnt - but would be nice to do this too ---- The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to see if your system is available in these configs please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation). -->

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended.

* `docker`
    * A generic configuration profile to be used with [Docker](https://docker.com/)
* `singularity`
    * A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
* `podman`
    * A generic configuration profile to be used with [Podman](https://podman.io/)
* `shifter`
    * A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
* `charliecloud`
    * A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
* `test`
    * A profile with a complete configuration for automated testing
    * Includes links to test data so needs no other parameters

### `-resume`

Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.
You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher requests (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

For example, if the nf-core/rnaseq pipeline is failing after multiple re-submissions of the `STAR_ALIGN` process due to an exit code of `137` this would indicate that there is an out of memory issue:

```console
[62/149eb0] NOTE: Process `RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137) -- Execution is retried (1)
Error executing process > 'RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)'

Caused by:
    Process `RNASEQ:ALIGN_STAR:STAR_ALIGN (WT_REP1)` terminated with an error exit status (137)

Command executed:
    STAR \
        --genomeDir star \
        --readFilesIn WT_REP1_trimmed.fq.gz  \
        --runThreadN 2 \
        --outFileNamePrefix WT_REP1. \
        <TRUNCATED>

Command exit status:
    137

Command output:
    (empty)

Command error:
    .command.sh: line 9:  30 Killed    STAR --genomeDir star --readFilesIn WT_REP1_trimmed.fq.gz --runThreadN 2 --outFileNamePrefix WT_REP1. <TRUNCATED>
Work dir:
    /home/pipelinetest/work/9d/172ca5881234073e8d76f2a19c88fb

Tip: you can replicate the issue by changing to the process work dir and entering the command `bash .command.run`
```

To bypass this error you would need to find exactly which resources are set by the `STAR_ALIGN` process. The quickest way is to search for `process STAR_ALIGN` in the [nf-core/rnaseq Github repo](https://github.com/nf-core/rnaseq/search?q=process+STAR_ALIGN). We have standardised the structure of Nextflow DSL2 pipelines such that all module files will be present in the `modules/` directory and so based on the search results the file we want is `modules/nf-core/modules/azimuth/main.nf`. If you click on the link to that file you will notice that there is a `label` directive at the top of the module that is set to [` label 'process_high_memory'`](https://github.com/wtsi-hgi/yascp/blob/73aed20e2c1b3ae693026704f2ccd3bb5d59deae/modules/nf-core/modules/azimuth/main.nf#L3). The [Nextflow `label`](https://www.nextflow.io/docs/latest/process.html#label) directive allows us to organise workflow processes in separate groups which can be referenced in a configuration file to select and configure subset of processes having similar computing requirements. The default values for the `process_high_memory` label are set in the pipeline's [`base.config`](https://github.com/wtsi-hgi/yascp/blob/73aed20e2c1b3ae693026704f2ccd3bb5d59deae/conf/base.conf#L200-L203) which in this case is defined as 50GB. Providing you haven't set any other standard nf-core parameters to __cap__ the [maximum resources](https://nf-co.re/usage/configuration#max-resources). Dependant on your setup in your institute you may want to overite these default parameters as per bellow: https://github.com/wtsi-hgi/yascp/blob/73aed20e2c1b3ae693026704f2ccd3bb5d59deae/conf/base.conf#L215-L218

```nextflow
process {
    withName: AZIMUTH{
        maxRetries    = 3
        errorStrategy = { task.attempt > 2 ? 'ignore' : 'retry' }
        memory = { check_max( 50.GB * task.attempt, 'memory' ) }
    }
}
```

### Tool-specific options

For the ultimate flexibility, we have implemented and are using Nextflow DSL2 modules in a way where it is possible for both developers and users to change tool-specific command-line arguments (e.g. providing an additional command-line argument to the `AZIMUTH` process) as well as publishing options (e.g. saving files produced by the `AZIMUTH` process that aren't saved by default by the pipeline). In the majority of instances, as a user you won't have to change the default options set by the pipeline developer(s), however, there may be edge cases where creating a simple custom config file can improve the behaviour of the pipeline.

As mentioned at the beginning of this section it may also be necessary for users to overwrite the options passed to modules to be able to customise specific aspects of the way in which a particular tool is executed by the pipeline. Given that all of the default module options are stored in the pipeline's `modules.config` as a [`params` variable](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/modules.config#L24-L25) it is also possible to overwrite any of these options via a custom config file.

Say for example we want to append an additional, non-mandatory parameter (i.e. `--outFilterMismatchNmax 16`) to the arguments passed to the `STAR_ALIGN` module. Firstly, we need to copy across the default `args` specified in the [`modules.config`](https://github.com/wtsi-hgi/yascp/blob/main/conf/modules.config) and create a custom config file that is a composite of the default `args` as well as the additional options you would like to provide. This is very important because Nextflow will overwrite the default value of `args` that you provide via the custom config.

<!-- WOULD BE NICE TO DO THIS - WE DONT HAVE IT CURRENTLY:
As you will see in the example below, we have:

* appended `--outFilterMismatchNmax 16` to the default `args` used by the module.
* changed the default `publish_dir` value to where the files will eventually be published in the main results directory.
* appended `'bam':''` to the default value of `publish_files` so that the BAM files generated by the process will also be saved in the top-level results directory for the module. Note: `'out':'log'` means any file/directory ending in `out` will now be saved in a separate directory called `my_star_directory/log/`.

```nextflow
params {
    modules {
        'star_align' {
            args          = "--quantMode TranscriptomeSAM --twopassMode Basic --outSAMtype BAM Unsorted --readFilesCommand zcat --runRNGseed 0 --outFilterMultimapNmax 20 --alignSJDBoverhangMin 1 --outSAMattributes NH HI AS NM MD --quantTranscriptomeBan Singleend --outFilterMismatchNmax 16"
            publish_dir   = "my_star_directory"
            publish_files = ['out':'log', 'tab':'log', 'bam':'']
        }
    }
}
``` -->

### Updating containers

The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. If for some reason you need to use a different version of a particular tool with the pipeline then you just need to identify the `process` name and override the Nextflow `container` definition for that process using the `withName` declaration. For example, in the [nf-core/viralrecon](https://nf-co.re/viralrecon) pipeline a tool called [Pangolin](https://github.com/cov-lineages/pangolin) has been used during the COVID-19 pandemic to assign lineages to SARS-CoV-2 genome sequenced samples. Given that the lineage assignments change quite frequently it doesn't make sense to re-release the nf-core/viralrecon everytime a new version of Pangolin has been released. However, you can override the default container used by the pipeline by creating a custom config file and passing it as a command-line argument via `-c custom.config`.

1. Check the default version used by the pipeline in the module file for [Pangolin](https://github.com/wtsi-hgi/yascp/blob/73aed20e2c1b3ae693026704f2ccd3bb5d59deae/modules/nf-core/modules/normalise_and_pca/main.nf#L15-L21)
2. Find the latest version of the Biocontainer available: [We currently dont have a default place to store these - need an s3 bucket]
3. Create the custom config accordingly:

    * For Docker:

        ```nextflow
        process {
            withName: PANGOLIN {
                container = 'quay.io/biocontainers/pangolin:3.0.5--pyhdfd78af_0'
            }
        }
        ```

    * For Singularity:

        ```nextflow
        process {
            withName: NORMALISE_AND_PCA {
                container = '/software/hgi/containers/wtsihgi_nf_scrna_qc_6bb6af5-2021-12-23-3270149cf265.sif'
            }
        }
        ```

    * For Conda:

        ```nextflow
        process {
            withName: PANGOLIN {
                conda = 'bioconda::pangolin=3.0.5'
            }
        }
        ```

> **NB:** If you wish to periodically update individual tool-specific results (e.g. Pangolin) generated by the pipeline then you must ensure to keep the `work/` directory otherwise the `-resume` ability of the pipeline will be compromised and it will restart from scratch.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```console
NXF_OPTS='-Xms1g -Xmx4g'
```
