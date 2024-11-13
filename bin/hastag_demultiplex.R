#!/usr/bin/env Rscript

library(Seurat)
library(tidyverse)

# Define base directory and specific suffixes for each data type
base_dir <- "."
cmo_suffix <- "__Multiplexing_Capture.tsv"
gex_suffix <- "__Gene_Expression"
adt_suffix <- "__Antibody_Capture"

# Find files or directories ending with the required suffixes
cmo_file <- list.files(base_dir, pattern = paste0(".*", cmo_suffix, "$"), full.names = TRUE)
gex_dir <- list.files(base_dir, pattern = paste0(".*", gex_suffix, "$"), full.names = TRUE, include.dirs = TRUE)
adt_dir <- list.files(base_dir, pattern = paste0(".*", adt_suffix, "$"), full.names = TRUE, include.dirs = TRUE)

# Check if files/directories exist
cmo_exists <- length(cmo_file) > 0
gex_exists <- length(gex_dir) > 0
adt_exists <- length(adt_dir) > 0

# Load Gene Expression data if it exists
if (gex_exists) {
    data_GEX <- Read10X(gex_dir)
    seurat <- CreateSeuratObject(counts = data_GEX, assay = "RNA")
} else {
    stop("Gene Expression data not found. Cannot proceed.")
}

# Load Antibody Capture data if it exists
if (adt_exists) {
    data_ADT <- Read10X(adt_dir)
    seurat[["ADT"]] <- CreateAssayObject(counts = data_ADT, assay = "ADT")
}

# Load CMO data if it exists, otherwise skip CMO-related processing
if (cmo_exists) {
    data_HTO <- read.csv(cmo_file, row.names = 1, sep = '\t')
    seurat[['CMO']] <- CreateAssayObject(counts = as(t(as.matrix(data_HTO)), "dgCMatrix"), assay = 'CMO')
    seurat <- NormalizeData(seurat, assay = 'CMO', normalization.method = "CLR", margin = 2)
    
    # Run MULTIseqDemux if CMO data is available
    seurat_object_demux <- MULTIseqDemux(
        object = seurat,
        assay = "CMO",
        quantile = 0.7,
        autoThresh = FALSE,
        maxiter = 5,
        qrange = seq(from = 0.1, to = 0.9, by = 0.05),
        verbose = TRUE
    )

    # Save results if demultiplexing was successful
    results <- seurat_object_demux@meta.data
    write.table(results, file = "hastag_demux_results.tsv", sep = "\t", row.names = TRUE, col.names = NA, quote = FALSE)
} else {
    message("CMO data not found. Skipping CMO-related processing.")
}