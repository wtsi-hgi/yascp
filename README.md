<p align="center">
  <img src="assets/images/YASCP_Logo.png" width="50%"/>
</p>

## Introduction
Usage [DOCUMENTATION](https://maxozo.github.io/yascp/)
<p align="center">
  <img src="https://github.com/wtsi-hgi/yascp/blob/main/assets/images/yascp_workflow.png" width="100%"/>
</p>

<p align="center">
  <img src="https://github.com/wtsi-hgi/yascp/blob/main/assets/illustrator_files/1x/Yascp_workflow-03.png" width="100%"/>
</p>

As indicated above you can run pipeline blocks independently:

<img width="100%" alt="Screenshot 2024-06-03 at 17 01 01" src="https://github.com/wtsi-hgi/yascp/assets/22347136/c724f731-42ab-4880-9666-eeb3384fd5e6">

**YASCP** (Yet Another Single Cell Pipeline) is a scalable and modular single-cell analysis pipeline designed for high-quality preprocessing, deconvolution, doublet detection, clustering, cell type assignment, and integration. The acronym moves from Y to A to symbolize the pursuit of knowledge — no Z implies there's always more to explore, refine, and improve.


YASCP supports:

- 10x Genomics Cell Ranger-style output directories 
- Raw BAM/BAI and MTX file inputs for nonstandard or intermediate-stage datasets 
- Hashtag multiplexing (HTO/CITE-seq)  
- CITE-seq protein expression quantification  
- scRNA-seq and scATAC-seq analysis modes  
- Modular execution: run the full pipeline or reuse individual steps independently
  
YASCP is not hardcoded for any specific tissue, platform, or cell type. Each module — from ambient RNA removal to clustering and annotation — can be configured independently, allowing you to tailor thresholds, skip irrelevant steps, or restart from any stage of the analysis.

### Flexibility for Custom Designs

YASCP is designed to accommodate complex experimental scenarios such as stimulation conditions, CITE-seq profiling, and multimodal assays. While it does not perform CRISPR guide assignment or perturbation modeling directly, it can preprocess CRISPR-based single-cell data upstream of specialized tools — providing normalized, QC-filtered .h5ad or .rds files ready for downstream integration.
Whether you’re analyzing a simple 10x run or building a customized, multi-step workflow, YASCP enables reproducible, flexible, and scalable analysis.

- **Condition-aware and hashtag-aware workflows**  
  You can split samples by hashtag or stimulation *before* QC, doublet detection, or annotation, either externally or by modifying the workflow schema.

- **Custom QC and filtering per sample/tag**  
  Per-donor or per-hashtag QC thresholds can be applied dynamically using modular blocks.

- **Optional or replaceable modules**  
  All major steps (e.g., CellBender, deconvolution, cell typing, doublet detection) can be skipped. You can use filtered CellRanger outputs directly.

- **Antibody and hashtag splitting**  
  CITE-seq protein and hashtag features can be split and analyzed separately.

- **Custom inputs at any stage**  
  You can re-enter the pipeline with intermediate `.h5ad` objects and resume downstream analysis.

- **Manual thresholds and logic injection**  
  Flexible logic and parameter overrides allow manual thresholding or customized QC/annotation rules per batch or condition.


Results will demultiplex individuals, robustly assess the assignments
![Screenshot 2024-06-03 at 12 56 44](https://github.com/wtsi-hgi/yascp/assets/22347136/5129c789-fbe9-41e8-8d28-5d286896f14a)


As well as assign celltypes, perform integrations, remove ambient RNA and produce publication ready plots
![Screenshot 2024-06-02 at 15 20 29](https://github.com/wtsi-hgi/yascp/assets/22347136/fe39d33a-97ec-44a1-9614-55f3585bde4d)

Developed by M.Ozols under the leadership of N.Soranzo and Human Genetics Informatics (HGI), this large-scale single-cell pipeline was originally crafted for the Cardinal project (profiling UKBB and GH participants) but is versatile enough for broad scRNA analysis applications. 

Input requires a tsv seperated file [(please read detailed documentation here)](https://github.com/wtsi-hgi/yascp/tree/yascp_docs) with paths and if running in an genotype  additional input is required to be provided in an input.nf file pointing to the vcf location. This pipeline is designed to be used any large scale single cell experiments.

The foundational ideas were inspired by earlier pipelines from Anderson lab but has been expanded, specifically those for [deconvolution](https://github.com/wtsi-hgi/nf_scrna_deconvolution.git), [cellbender](https://github.com/wtsi-hgi/nf_cellbender), and [quality control and clustering](https://github.com/wtsi-hgi/nf_qc_cluster/tree/main). This ensures a robust integration of proven methodologies tailored to meet the demands of expansive single-cell data analysis.


<!-- TODO nf-core: Add full-sized test dataset and amend the paragraph below if applicable -->
<!-- On release, automated continuous integration tests run the pipeline on a full-sized dataset on the AWS cloud infrastructure. This ensures that the pipeline runs on AWS, has sensible resource allocation defaults set to run on real-world datasets, and permits the persistent storage of results to benchmark between pipeline releases and other analysis sources. The results obtained from the full-sized test can be viewed on the [nf-core website](https://nf-co.re/yascp/results). -->

## Quick Start

1. Install [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation) (`>=21.04.0`)
2. Install any of [`Docker`](https://docs.docker.com/engine/installation/), [`Singularity`](https://www.sylabs.io/guides/3.0/user-guide/)for full pipeline reproducibility.
3. Download/clone the pipeline and test it on a minimal dataset with a single command:

    ```console
    git clone https://github.com/wtsi-hgi/yascp.git
    nextflow run /path/to/colned/yascp -profile test,<docker/singularity,institute>
    ```
## Run on Your Own Data

  ### 1. Prepare an `input.tsv` file

  You can run YASCP in two ways depending on the structure of your input data.

  #### Option 1: Minimal (standard 10x format)

  Use this format if you have 10x Genomics-style output folders (e.g., from Cell Ranger):

  | experiment_id | n_pooled | donor_vcf_ids        | data_path_10x_format        |
  |---------------|----------|-----------------------|-----------------------------|
  | SampleA       | 1        | ""                    | /data/project1/10x_output1/ |
  | SampleB       | 2        | D001,D002             | /data/project1/10x_output2/ |

  > `donor_vcf_ids` should be comma-separated if multiple donors are pooled. Leave empty ("") if not applicable.

  ---

  #### Option 2: Custom paths (nonstandard structure or preprocessed inputs)

  Use this format if you want to supply filtered/unfiltered MTX files or BAM/BAI directly:

  | experiment_id | n_pooled | donor_vcf_ids  | data_path_10x_format | filtered_mtx             | filtered_barcodes         | filtered_features         | unfiltered_mtx           | unfiltered_barcodes       | unfiltered_features       | bam                     | bai                     |
  |---------------|----------|----------------|-----------------------|---------------------------|----------------------------|----------------------------|---------------------------|----------------------------|----------------------------|--------------------------|--------------------------|
  | SampleC       | 4        | D101,D102,D103 | ""                    | /data/SampleC/filtered.mtx.gz | /data/SampleC/filtered.barcodes.tsv.gz | /data/SampleC/filtered.features.tsv.gz | /data/SampleC/unfiltered.mtx.gz | /data/SampleC/unfiltered.barcodes.tsv.gz | /data/SampleC/unfiltered.features.tsv.gz | /data/SampleC/data.bam | /data/SampleC/data.bam.bai |
  | SampleD       | 3        | D201,D202,D203 | ""                    | ...                       | ...                        | ...                        | ...                       | ...                        | ...                        | ...                      | ...                      |

  - If `data_path_10x_format` is provided, it takes precedence.
  - If empty (`""`), the pipeline will fall back to the provided `filtered_*`, `unfiltered_*`, or `bam`/`bai` file paths.
  - You must provide either a valid 10x directory or a complete alternative set.

  ---

  ### 2. Run the pipeline

  ```bash
  git clone https://github.com/wtsi-hgi/yascp.git
  cd yascp

  nextflow run ./main.nf \
    -profile <test,docker,singularity,institute> \
    --input_data_table path/to/input.tsv
  ```

## Pipeline summary
Pipeline has a modular design ensuring that the bits and piecies can be run independently according to project needs. Overall pipeline is focussed arounf main steps:
1. Cellbender
2. CellSNP
3. Vireo
4. Souporcell
5. Celltypist
6. Azimuth
7. BBKNN
8. Harmony
9. Scrublet, DoubletDecon, DoubletFinder, SCDS, scDblFinder, DoubletDetection
10. Sccaf
11. Lisi
12. Isolation Forest
13. Hard filters
14. Genotype deconvolution and GT match against multiple panels.
15. Citeseq DSB normalisations, 
16. Cell genotype concordance Calculations




## Documentation: Prepearing your own data and interpreting the results

To understand how to prepeare your own data and how to interpret the results please refear to [documents HERE](https://github.com/wtsi-hgi/yascp/tree/yascp_docs)


## Credits

Yascp was originally written by Matiss Ozols with major contributions from Leland Taylor, Guillaume Noell, Hannes Ponstingl, Vladimir Ovchinnikov,  Vivek Iyer,  Henry Taylor, Tobi Alegbe, Monika Krzak, Alessandro Raveane, Carl Anderson, Anna Lorenc, Haerin Jang, Niek de Klein, Stephen Watt, Nicole Soranzo, Oliver Stegle.

The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across multiple compute infrastructures in a very portable manner. It uses Docker/Singularity containers making installation trivial and results highly reproducible. The [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) implementation of this pipeline uses one container per process which makes it much easier to maintain and update software dependencies. Where possible, these processes have been submitted to and installed from [nf-core/modules](https://github.com/nf-core/modules) in order to make them available to all nf-core pipelines, and to everyone within the Nextflow community!

## Contributions and Support

We wellcome all contributions. If you would like to contribute to this pipeline, please create a fork and then create a pull request, and inform Matiss (mo11@sanger.ac.uk) re the changes made and additions added.

## Citations


If you use  YASCP for your analysis, please cite it using the following doi: [10.5281/zenodo.15600242 ](https://doi.org/10.5281/zenodo.15600242 )
> Ozols, M. et al. YASCP: GitHub. https://github.com/wtsi-hgi/yascp [![DOI](https://zenodo.org/badge/423515781.svg)](https://doi.org/10.5281/zenodo.15600242)
We are also working on publishing this pipeline. 

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

We have used nf-cores template to develop this pipeline. You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
