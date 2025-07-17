#!/usr/bin/env Rscript

## expects input *.h5 file as argument
## writes plots to a single PDF file Rplots.pdf
library(Azimuth)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(hdf5r)
library(ggplot2)
library(tools)
library(Signac)
library(Seurat)
library(GenomicRanges)
options(future.globals.maxSize = 2000 * 1024^2)
args =list()
# inputfile.h5ad='./Pool1.h5ad'
# REFERENCE_DIR='/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/jaguar_yascp/nieks_pipeline/yascp_run/ref_kidney'
# levels='annotation.l3,annotation.l2,annotation.l1'

args = commandArgs(trailingOnly=TRUE)


inputfile.h5ad = args[1]
REFERENCE_DIR = args[2]

if (REFERENCE_DIR=='PBMC'){
  REFERENCE_DIR <- "/opt/PBMC_reference/"
}

levels = args[3]
prefix = args[4]
# levels = "celltype.l2,celltype.l1,celltype.l3"
levels = unlist(x = strsplit(x = levels, split = ',', fixed = TRUE))
query_matrix <- Seurat::Read10X(inputfile.h5ad)
query <- CreateSeuratObject(counts = query_matrix)
query$nCount_RNA <- Matrix::colSums(GetAssayData(query, assay = "RNA", layer = "counts"))


peak_annot <- read.table(paste0(inputfile.h5ad,'/',"peaks.bed"), sep = "\t", header = FALSE)
colnames(peak_annot) <- c("chr", "start", "end")
peaks <- makeGRangesFromDataFrame(peak_annot)
# Create Chromatin Assay
# Fully normalize the path and force bgzip check manually
fragment_path <- normalizePath(file.path(inputfile.h5ad, "fragments.tsv.gz"))
h5_path <- normalizePath(file.path(inputfile.h5ad, "filtered_peak_bc_matrix.h5"))
print(fragment_path)
# Optional: test tabix index validity before proceeding
if (!file.exists(fragment_path)) stop("fragments.tsv.gz not found.")
if (!file.exists(paste0(fragment_path, ".tbi"))) stop("fragments.tsv.gz.tbi index not found.")

# # Create Fragment object from path
# frag_obj <- CreateFragmentObject(
#   path = fragment_path,
#   cells = colnames(query_matrix)  # important: restrict to observed barcodes
# )

# # Now make the chromatin assay
# chrom_assay <- CreateChromatinAssay(
#   counts = query_matrix,
#   sep = c(":", "-"),
#   fragments = frag_obj,
#   ranges = peaks
# )

# # Create Seurat object
# atac <- CreateSeuratObject(counts = chrom_assay, assay = "peaks")

# The RunAzimuth function can take a path or Seurat object as input
atac <- RunAzimuth(query = h5_path, query.modality = "ATAC",
    reference = '/lustre/scratch127/humgen/teams_v2/hgi/mo11/tmp_projects/cardinal/cardinal_atac/develop/v1/work/38/170f2e5434d4bd58aae30244c3a025/45c38c611dca591a842a70ce621a7081__Gene_Expression', fragment.path = fragment_path)

metadata <- atac_res@meta.data
write.csv(metadata, file = "atac_azimuth_metadata.csv", row.names = TRUE)
