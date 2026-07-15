# nf-core/yascp: Output

This document describes the output folder structure produced by a full YASCP pipeline run. By default, results are written to `results/` (configurable via `--outdir`). Not all directories will be present in every run — the pipeline produces outputs only for the steps you enable.

## Directory Overview
----
```
results/
├── preprocessing/              # Data preparation
├── deconvolution/              # Donor deconvolution & genotype matching
├── doublet_detection/          # Doublet identification
├── celltype_assignment/        # Cell type annotation
├── handover/                   # Final merged datasets & summary reports
├── citeseq/                    # CITE-seq / multi-modal integration
├── clustering_and_integration  # Clustering and integration
├── plots/                      # Plots
├── pipeline_info/              # Nextflow execution reports
└── yascp_software_versions.yml # Software versions used in the run
```
## 1. preprocessing
```
preprocessing/
├── cellbender/                           # CellBender corrected matrices & QC plots per sample
├── resources/                            # Downloaded reference resources
├── data_modalities_split/{mode}/{pool}/  # GEX/ADT/HTO split
└── subset_genotypes/                     # Subsetted VCFs per pool
```
## 2. deconvolution
```
deconvolution/
├── cellsnp/cellsnp_{pool}/                    # CellSNP pileup outputs per pool
├── mpileup/                                   # BAM pileup (BAM-derived SNP panels)
├── vireo/
│   ├── vireo_raw/{pool}/                      # Raw Vireo results
│   ├── vireo_processed/{pool}/                # Post-processed with corrected donor IDs
│   └── vireo_subsampling_cellsnp/{pool}/      # Subsampling robustness
├── split_donor_h5ad/{pool}/                   # Per-donor h5ad files
├── concordances/                              # Genotype concordance
├── gtmatch/                                   # Genotype matching & donor assignments
├── infered_genotypes/{pool}/                  # Inferred VCFs from Freebayes
└── hastag_demultiplex/{pool}/                 # HTO demultiplexing
```
## 3. doublet_detection
```
doublet_detection/
├── scrublet/                              # Scrublet results
├── multiplet.method=doubletdetection/     # DoubletDetection results
├── DoubletDecon/                          # DoubletDecon results
├── DoubletFinder/                         # DoubletFinder results
├── scDblFinder/                           # scDblFinder results
├── SCDS/                                  # SCDS results
├── doublet_results_combined/              # Combined results per pool
└── droplet_type_distribution/             # Distribution plots
```
## 4. celltype_assignment
```
celltype_assignment/
├── azimuth/{refset_name}/`         #Azimuth predictions per reference dataset
├── celltypist/{model}/{pool}/      #CellTypist predictions per model and sample
├── keras_celltype/{pool}/          #Keras deep-learning cell type predictions
├── scpred/                         #scPred cell type predictions
├── All_Celltype_Assignments.tsv    #File combining all cell type annotations across methods
├── donor_celltype_report.tsv       #Cell type summary per donor
└── tranche_celltype_report.tsv     #Cell type summary per tranche/batch
```
## 5. handover
```
handover/
├── merged_h5ad/`                   #Merged AnnData (h5ad) files at various pipeline stages
├── Donor_Quantification/{sample}/` #Per-donor BAM files and quantification
├── Donor_Quantification_summary    #Summary tables
├── Summary_plots                   #Summary plots
└── UMAPs                           #UMAP plots
```
## 6. citeseq
```
citeseq/
├── `DSB/{sample_name}/`   #DSB-normalised protein expression per sample
└── `all_data_integrated/` #Seurat WNN integrated multi-modal results, including VDJ integration if available
```
## 7. clustering_and_integration
```
clustering_and_integration/
├── plots                                  # plots
└── normalize=total_count.${param_details} # PCA results
```
