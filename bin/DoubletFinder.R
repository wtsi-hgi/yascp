#!/usr/bin/env Rscript

suppressMessages(suppressWarnings(library(argparse)))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-o", "--out", required = TRUE, help="The output directory where results will be saved")
parser$add_argument("-s", "--seurat_object", required = TRUE, type = "character", help = "A QC, normalized seurat object with classifications/clusters as Idents() saved as an rds object.")
parser$add_argument("-c", "--sct", required = TRUE, type = "logical", help = "Whether sctransform was used for normalization.")
# parser$add_argument("-d", "--doublet_number", required = TRUE, type = "integer", help = "Number of expected doublets based on droplets captured.")
parser$add_argument("-mr", "--expected_multiplet_rate", required = TRUE, type = "double", help = "Number of expected doublets rate based on droplets captured.")
parser$add_argument("-p", "--PCs", required = FALSE, default = 10, type = "integer", help = "Number of PCs to use for \'doubletFinder_v3\' function.")

parser$add_argument("-n", "--pN", required = FALSE, default = 0.25, type = "double", help = "Number of doublets to simulate as a proportion of the pool size.")

# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(DoubletFinder)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tidyr)))
suppressMessages(suppressWarnings(library(tidyverse)))

## make sure the directory exists ###
dir.create(args$out, recursive = TRUE)

## Add max future globals size for large pools

library(viridis)
library(SeuratDisk)
library(future)
options(future.globals.maxSize= 1020971520000)
if (future::supportsMulticore()) {
  future::plan(future::multicore)
} else {
  future::plan(future::multisession)
}

## make sure the directory exists ###
dir.create(args$out, recursive = TRUE, showWarnings = FALSE)

## Read in Data ##
# seurat <- readRDS(args$seurat_object)

Convert(
  args$seurat_object,
  dest = paste('tmp',"h5seurat",sep='.'),
  assay = "RNA",
  overwrite = TRUE,
  verbose = TRUE,
)
seurat <- LoadH5Seurat(paste('tmp',"h5seurat",sep='.'),assays = "RNA")
seurat <- NormalizeData(seurat)
print('Normalised')
all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)
print('Scaled')
seurat <- FindVariableFeatures(object = seurat)
# Check the number of variable features
num_samples <- ncol(seurat)
if (num_samples < 50) {
    print("No sufficinet number of cells to detect doublets.")
    quit(status = 0) 
}
# Determine the number of PCs to use
npcs_to_use <- ifelse(num_samples > 50, 50, num_samples-1)

# Run PCA with the determined number of PCs
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat), npcs = npcs_to_use)
print('PCA performed')
seurat <- FindNeighbors(seurat, dims = 1:10)
print('Neighbors found')
seurat <- FindClusters(seurat, resolution = 0.5)
print(seurat[["pca"]], dims = 1:5, nfeatures = 5)


### Read in the data
# seurat <- readRDS(args$seurat_object)


## pK Identification (no ground-truth) ---------------------------------------------------------------------------------------
sweep.res.list <- paramSweep(seurat, PCs = 1:10, sct = args$sct)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)
plot <- ggplot(bcmvn, aes(pK, BCmetric)) +
    geom_point()
ggsave(plot, filename = paste0(args$out,"/pKvBCmetric.png"))

## Homotypic Doublet Proportion Estimate -------------------------------------------------------------------------------------
annotations <- Idents(seurat)
homotypic.prop <- modelHomotypic(annotations)
# nExp_poi <- args$doublet_number
nExp_poi = round(args$expected_multiplet_rate*nrow(seurat@meta.data)) 
print(paste0("Expected number of doublets: ", nExp_poi))
nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))
## Run DoubletFinder with varying classification stringencies ----------------------------------------------------------------
seurat <- doubletFinder(seurat, PCs = 1:args$PCs, pN = args$pN, pK = as.numeric(as.character(bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))])), nExp = nExp_poi.adj, reuse.pANN = FALSE, sct = args$sct)
doublets <- as.data.frame(cbind(colnames(seurat), seurat@meta.data[,grepl(paste0("pANN_0.25_",as.numeric(as.character(bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))]))), colnames(seurat@meta.data))], seurat@meta.data[,grepl(paste0("DF.classifications_0.25_",as.numeric(as.character(bcmvn$pK[which(bcmvn$BCmetric == max(bcmvn$BCmetric))]))), colnames(seurat@meta.data))]))
colnames(doublets) <-  c("Barcode","DoubletFinder_score","DoubletFinder_DropletType")
doublets$DoubletFinder_DropletType <- gsub("Singlet","singlet",doublets$DoubletFinder_DropletType) %>% gsub("Doublet","doublet",.)

write_delim(doublets, file = paste0(args$out,"/DoubletFinder_doublets_singlets.tsv"), delim = "\t")

### Calculate number of doublets and singlets ###
summary <- as.data.frame(table(doublets$DoubletFinder_DropletType))
colnames(summary) <- c("Classification", "Droplet N")
write_delim(summary, paste0(args$out,"/DoubletFinder_doublet_summary.tsv"), "\t")