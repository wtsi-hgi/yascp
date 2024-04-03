# nf-core/yascp: Output

## Introduction

This document describes the output produced by the pipeline.

<!-- TODO nf-core: Write this documentation describing your workflow's output -->

# Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps:
The overall results folder will look simmillar to this:

![Screenshot 2024-04-02 at 15 08 13](https://github.com/wtsi-hgi/yascp/assets/22347136/12cc3575-8772-43ee-b64d-bb396e10ba82)

Where we have outputs from different steps of pipeline:
* [cellsnp](#cellsnp)
* [celltype identification](#celltype-identification)
* citeseq data processing
* [clustering and integration](#integration-and-clustering)
* [sample deconvolution](#vireo)
* [doublet detection](#doublet-detection)
* [genotype match](#vireo) to determine sample matches
* [handover](#handover) folder where summary statistics, plots and final qcd and annotated h5ads per donor are stored at.
* infered genotypes - output from vireo that has generated vcf files for each of the deconvoluted donors in pool.
* merged_h5ads - different preprocessing step merged h5ads (these allow to start the pipeline again in a clustering only mode)
* [nf-preprocessing](#ambient-rna-removal) - contains cellbender results
* pipeline info - statistics of the pipeline run.
* plots - some quality control plots.
* recourses - reference genome used in data processing.
* UMAPS - summary plot UMAPS - for a quick look. 

Each of these steps and the outputs produced are decribed more in detail bellow:

## Alignment step
#### [Cellranger](#Cellranger) - Curently users have to run Cellranger upstream of pipeline - we suggest to use the [no-cores pipeline](https://nf-co.re/scrnaseq/2.5.1) - https://nf-co.re/scrnaseq/2.5.1
## Ambient RNA removal
#### [Ambient RNA Removal using Cellbender](#Cellbender) - Reads the Cellranger outputs and removes the ambient RNA using [Cellbender](https://github.com/broadinstitute/CellBender)

<details markdown="1">
<summary>Output file structure ( nf-preprocessing/cellbender ):</summary>

*   Here we have multiple different plots and output files, however the most important ones are the matrix and h5ad files after the ambient rna removal: such as cellbenderFPR_0pt1filtered_10x_mtx/ cellbender_FPR_0.1_filtered.h5
    * ![Cellbender module output structure](../assets/images/cellbender_output_structure.png)
</details>

<details markdown="1">
<summary>Cellbender output plots:</summary>

*   Cellbender output plots:
    * ![Cellbender UMAP plot](../assets/images/cb_umap.png)
</details>

## Genotype processing and Donor deconvolutions (if more than 1 donor is in the pool) and Multiplet/Unassigned cell removal
####  [Genotype processing](#Genotype_processing) - If users provide the genotypes this step slices and dices the genotypes to prepeare these for the CellSNP/Vireo deconvolutions and GT matches
#### [Donor Deconvolution using CellSnp/Vireo](#CellSnp/Vireo) - We run cellsnp and vireo to deconvolute donors if the input file has indicated that there are more than 1 donors in the pool.

#### Cellsnp
<details markdown="1">
<summary>Cellsnp Output files:</summary>

* Cellsnp profiles each of the droplets for the variants in them, which is later utilised by vireo to assign the particular cell to the donor cluster:
    * ![Cellsnp output structure](../assets/images/cellsnp.png)
</details>

<details markdown="1">
<summary>Vireo Output files:</summary>

#### Vireo
* Vireo takes the cellsnp variant pileups and assigns donors the particular cell to the donor cluster:
    * ![Vireo output structure](../assets/images/Vireo_outputs.png)
</details>

#### Doublet Detection
![Screenshot 2024-04-02 at 15 43 16](https://github.com/wtsi-hgi/yascp/assets/22347136/781ce3b7-ea5e-4fe4-9ca3-d16e8b47123e)
<details markdown="1">
<summary>Scrublet Output files:</summary>

* By default we always run Scrublet - if we have no donors pooled in the run (i.e if we have only 1 donor), then the doublets will be removed by scrublet instead of vireo:
    * ![Scrublet output structure](../assets/images/Scrublet.png)
</details>

<details markdown="1">
<summary>DoubletDecon Output files:</summary>
    
* DoubletDecon output files contain barcode and label of whether its a singlet or a doublet:
    * ![Screenshot 2024-04-02 at 15 51 26](https://github.com/wtsi-hgi/yascp/assets/22347136/603d27e1-42e3-4be7-bbfd-ebb3412b3ec4)

</details>

<details markdown="1">
<summary>doubletdetection Output files:</summary>
    
* doubletdetection output files contain barcode and label of whether its a singlet or a doublet:
    * ![Screenshot 2024-04-02 at 15 59 15](https://github.com/wtsi-hgi/yascp/assets/22347136/c798d675-c96d-4137-92c2-6fa9340437c5)

</details>

<details markdown="1">
<summary>DoubletFinder Output files:</summary>
    
* DoubletFinder output files contain barcode and label of whether its a singlet or a doublet:
    * ![Screenshot 2024-04-02 at 16 00 47](https://github.com/wtsi-hgi/yascp/assets/22347136/4cdd8ba2-5d16-4c9b-a64e-aa9423514208)


</details>

<details markdown="1">
<summary>scDblFinder Output files:</summary>
    
* scDblFinder output files contain barcode and label of whether its a singlet or a doublet:
    * ![Screenshot 2024-04-02 at 16 01 49](https://github.com/wtsi-hgi/yascp/assets/22347136/69f7b19f-3b22-46bb-aafe-403ca3c399ae)

</details>


<details markdown="1">
<summary>SCDS Output files:</summary>
    
* SCDS output files contain barcode and label of whether its a singlet or a doublet:
    * ![Screenshot 2024-04-02 at 16 02 23](https://github.com/wtsi-hgi/yascp/assets/22347136/b2b8ca81-449b-4a94-a1a9-2bec592e74f4)

</details>

#### [Donor Deconvolution using Souporcell](#Souporcell) - Souporcell option both removes the ambioent RNA and deconvolutes the donors [currently however this option is broken and will be fixed soon]

#### [GT match](#GT_match) - This step utilises the prepeared genotypes and the infered genotypes by Vireo and picks out the donor that corresponds to the right reads. 

<details markdown="1">
<summary>GT input files:</summary>

* Users can provide multipple different cohort VCFs and that are split per chromosomes or one big vcf/bcf file.:
    * ![GT input structure](../assets/images/gt_input_vcfs.png)
</details>

<details markdown="1">
<summary>GT match results structure:</summary>

* GT match produces multiple metrics that assesses whether donor is the one we expect and what is the relatedness within pool.
    * ![GT input structure](../assets/images/gt_structure.png)

* Results indicate which donor from Vireo deconvolutions is which:
    * ![GT input structure](../assets/images/GT_Results.png)
</details>


## Celltype identification
#### [Azimuth](#Azimuth) - Uses Azimuth PBMC l2 reference (pipeline will be adjusted later to be more general for other tissue types) to assign the celltypes. Downstream it maps the l2 to l1 and l3 as per https://github.com/wtsi-hgi/yascp/blob/main/assets/azimuth/Azimuth_Mappings.txt 

<details markdown="1">
<summary>Azimuth Output files:</summary>

* By default we run azimuth l2 celltype assignment:

    * ![Scrublet output structure](../assets/images/Azimuth.png)
</details>

#### [Celltypist](#Celltypist) - Performs cellype assignment using celltypist Imule Low and Imune High profiles (this will be adjusted to use more references)

<details markdown="1">
<summary>Celltypist Output files:</summary>

* By default we run Imune High, Imune Low and Imune PBMC reference celltype assignment:

    * ![Celltypist output structure](../assets/images/Celltypist.png)
</details>

<details markdown="1">
<summary>Combined celltypes file:</summary>


#### [Keras celltype transfer](#Keras) - This is utilising pretrained reference panels for celltype assignment - curently only works in Sanger.

#### Combined File - A combined Celltypes file is produced by pipeline where all different references are combined in one spreadsheet.:

    * ![Celltypist output structure](../assets/images/Combined_celltypes.png)
</details>

## Donor and Cell QC
### We perform different types of QC, Adaptive Isolation Forests, Adaptive Isolation Forests per celltype, Hard Filters tresholds.
<details markdown="1">
<summary>Data QC output folder structure:</summary>

*   QC output Folder structure:
    * ![Clustering BBKNN structure](../assets/images/QC_structure.png)
</details>

####  [Isolation Forest](#Isolation_Forest)
<details markdown="1">
<summary>We parfor Isolation forests in different resolutions - All data together, Per Celltype adaptive qc:</summary>

*   All together Isolation Forests:
    * ![Adaptive structure](../assets/images/adaptive_qc_alltogether.png)
*   Per Celltype Isolation Forests:
    * ![Adaptive structure](../assets/images/adaptive_qc_celltype.png)
</details>

####  [Hard filters](#Hard_filters): We also perform hard filters if user has specified that this is something thats required.

## Integration and clustering
### By default multiple different clustering resolutions will be run for both BBKNN and Harmony resulting in a subfolder structure. Pipeline automatically estimates the best number of PCs to use for clustering using knee and elbow plots that can be found in plots section.
<details markdown="1">
<summary>Output file structure ( clustering ):</summary>

*   Clustering combines all different integration methodologies utilised and in addition different plots in a structure represented in this layout:
    * ![Clustering module output structure](../assets/images/Clustering.png)
</details>


#### [BBKNN](#BBKNN) - 

<details markdown="1">
<summary>BBKNN file structure ( clustering ):</summary>

*   BBKNN is performed with different clustering resolutions and each of the clusters assesed ussing sccaf:
    * ![Clustering BBKNN structure](../assets/images/bbknn_structure.png)
</details>

<details markdown="1">
<summary>BBKNN sample UMAPS Coloured:</summary>

*   Resolution 0.1: BBKNN is performed with different clustering resolutions and each of the clusters assesed ussing sccaf:
    * ![Clustering BBKNN structure](../assets/images/bbknn_cluster1.png)
*   Resolution 5: BBKNN is performed with different clustering resolutions and each of the clusters assesed ussing sccaf:
    * ![Clustering BBKNN structure](../assets/images/bbknn_cluster2.png)
*   Mitochondial transcripts: Coloured UMAP: We also color each of the bespoke clusters with different metrics:
    * ![Clustering BBKNN structure](../assets/images/BBKNN_mito.png)
</details>

#### [Harmony](#Harmony) - 

<details markdown="1">
<summary>Harmony file structure ( clustering ):</summary>

*   Harmony is performed with different clustering resolutions and each of the clusters assesed ussing sccaf:
    * ![Clustering Harmony structure](../assets/images/harmony_structure.png)
</details>
<details markdown="1">
<summary>Harmony sample UMAPS Coloured:</summary>

*   Resolution 0.1: Harmony is performed with different clustering resolutions and each of the clusters assesed ussing sccaf:
    * ![Clustering Harmony structure](../assets/images/harmony_cluster1.png)
*   Resolution 5: Harmony is performed with different clustering resolutions and each of the clusters assesed ussing sccaf:
    * ![Clustering Harmony structure](../assets/images/harmony_cluster2.png)
*   Mitochondial transcripts: Coloured UMAP: We also color each of the bespoke clusters with different metrics:
    * ![Clustering Harmony structure](../assets/images/harmony_mito.png)
</details>

<details markdown="1">
<summary>Harmony cluster evaluations and cluster markers:</summary>

*   Histograms: Multiple useful prolts are produced to look at the clusterings:
    * ![Harmony Histogram](../assets/images/harmony_histo.png)
*   Dotplots: Multiple useful prolts are produced to look at the clusterings:
    * ![Harmony Dotplot](../assets/images/harmony_dotplot)
</details>

#### [PCA](#BBKNN) 

<details markdown="1">
<summary>PCA file structure ( clustering ):</summary>

*   PCA is performed on the integrated data:
    * ![Clustering Harmony structure](../assets/images/pca_plots.png)
</details>
<details markdown="1">
<summary>PCA file structure ( clustering ):</summary>

*   Gene Loadings for each of the PCA is evaluated:
    * ![PCA gene loadings](../assets/images/PCA_Gene_loadings.png)
</details>



## Cluster assesments
#### [Sccaf](#Sccaf) We perform Sccaf to asses the clustering accuracies, these are useful metrics in picking the best resolution for clustrering.
<details markdown="1">
<summary>Sccaf file structure ( clustering ):</summary>

*   As described above clustering is assesed using scaff: directory structure:
    * ![Clustering BBKNN structure](../assets/images/sccaf_structure.png)

*   Precission recall curves:
    * ![Clustering BBKNN structure](../assets/images/scaff_prec_recal.png)

*   ROC:
    * ![Clustering BBKNN structure](../assets/images/sccaff_roc.png)

*   Accuracy:
    * ![Clustering BBKNN structure](../assets/images/sccaff_accuracy.png)

</details>

#### [Lisi](#Lisi) We also have a capability in running LISI cluster assesments, however curently this option does not run by default as it is memory demanding and requires some further optimisations

## [Handover](#handover): Summary Statistics, Per Donor h5ad files, Summary Plots
    * ![Screenshot 2024-04-03 at 16 44 35](https://github.com/wtsi-hgi/yascp/assets/22347136/64bd3ca8-cb10-48bb-8334-f12482e4ebfd)
In this folder we can see 3 different folders:
* Donor_Quantification - where we can see the Cellranger filtered, Cellranger raw, Cellbender filtered files that are used to produce the filal per donor h5ad files and the metadata features in the per donor tsv files
    * ![Screenshot 2024-04-03 at 16 47 28](https://github.com/wtsi-hgi/yascp/assets/22347136/8524243f-4bf1-4713-9076-0e1d3fcb99e1)

* Donor_Quantification_summary folder where we have summary statistics per donor and summary statistics per tranche (collection of all pools that were run in this run).
    * ![Screenshot 2024-04-03 at 16 49 59](https://github.com/wtsi-hgi/yascp/assets/22347136/a6107709-f83e-45e1-9008-1cdde1510c67)

* Summary _plots contains the most important plots per each of the steps for a quick inversigations of the performance of the scRNA runs and the performance of the analysis.
    * ![Screenshot 2024-04-03 at 16 51 35](https://github.com/wtsi-hgi/yascp/assets/22347136/7c63e2c0-6251-4a7d-8e14-be434c0e017b)



[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.



