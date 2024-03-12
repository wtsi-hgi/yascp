#!/usr/bin/env Rscript

# This script makes a seurat object for each sample (that has CITE measurements) 
# containing 3 assays: RNA, CITE, CITE_bgRemoved
# For CITE_bgRemoved assay, background (and isotype levels) are removed using the
# dsb package. 
# .libPaths('/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/jaguar_yascp/nieks_pipeline/yascp_run/work/d4/3825bd8202a0dfc28c0fffdb312c5e/libs')
# install.packages("ggpubr")
# install.packages("viridis")
#### libraries ####
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(ggpubr)
library(dsb)
library('Seurat')
library(viridis)
cat('Libraries are loaded')
library(SeuratDisk)
library(future)
options(future.globals.maxSize= 1020971520000)
if (future::supportsMulticore()) {
  future::plan(future::multicore)
} else {
  future::plan(future::multisession)
}

#####
args = commandArgs(trailingOnly=TRUE)
#### set up directories, colors paths ####
data_dir <- getwd()
outdir <- getwd()




#####

#### Get location of files with filtered cellranger .h5 files #### 
# The reason that recursive search of files is used is because different
# version of cellranger have different file paths and file output names
# This has been tested for Cellranger version 6 and 7, might not work with
# older or newer version

# get the names of the filtered feature files (from 'cellranger count')
# filtered_feature_file = cellranger_filepath = args[2]
#    
# filtered_cellranger = '/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/jaguar_yascp/nieks_pipeline/fetch/results_old/cellranger_data/cellranger700_multi_bc45a1c2fe2a3fbbcde46cf984cf42e2/per_sample_outs/cellranger700_multi_bc45a1c2fe2a3fbbcde46cf984cf42e2/count/sample_filtered_feature_bc_matrix.h5'
# sample_name <- 'cellranger700_multi_d00ddfc9ace09586f1ac108fca92db61'
# cellranger_rawfile_path <- 'raw_feature_bc_matrix'
# filtered_cellranger = 'filtered_feature_bc_matrix'
# file_with_qc_applied = 'cellranger700_multi_d00ddfc9ace09586f1ac108fca92db61___sample_QCd_adata.h5ad'
sample_name <- args[1]
cellranger_rawfile_path <- args[2]
filtered_cellranger = args[3]
file_with_qc_applied = args[4]
# cellranger_rawfile_path = '/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/jaguar_yascp/nieks_pipeline/fetch/results_old/cellranger_data/cellranger700_multi_bc45a1c2fe2a3fbbcde46cf984cf42e2/multi/count/raw_feature_bc_matrix.h5'

# sample_name <-'cellranger700_multi_bc45a1c2fe2a3fbbcde46cf984cf42e2'

# file_with_qc_applied <-'/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/jaguar_yascp/nieks_pipeline/yascp_run/work/a0/dd1ce16989fca8ff81e0d759589a1f/donor_level_anndata/cellranger700_multi_bc45a1c2fe2a3fbbcde46cf984cf42e2___sample_QCd_adata.h5ad'



# Create all output/figure directories
dir.create(outdir,showWarnings = F)
# dir.create(paste0(figdir,'/ADT_background'), showWarnings = F)
tmp_rds_dir <- paste0(outdir,'/tmp_rds_files/')
dir.create(tmp_rds_dir,showWarnings = F)
# For saving intermediary files
tmp_rds_dir <- paste0(outdir,'/tmp_rds_files__',sample_name,'/')
dir.create(tmp_rds_dir,showWarnings = F)
dir.create(paste0(tmp_rds_dir,'/1.CITE/'), showWarnings = F)

cite_data_dir <- paste0(outdir,'/CITE__',sample_name,'/')
dir.create(cite_data_dir,showWarnings = F)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))




Convert(
  file_with_qc_applied,
  dest = paste('tmp',"h5seurat",sep='.'),
  assay = "RNA",
  overwrite = TRUE,
  verbose = TRUE,
)
# for 'cellranger multi' name is different

# sample_feature_files <- list.files(path=data_dir,
#                                    pattern = "sample_filtered_feature_bc_matrix.h5$", 
#                                    recursive = TRUE, full.names=T)
# print(paste0(length(c(filtered_feature_files,sample_feature_files)),' h5 files found'))
#####
#####

#### Get samplesheets for all batches and etact sample IDs ####
# Probably does not need 
# samplesheet <- data.frame()
# sampleID_files <- list.files(path = paste0(Sys.glob(paste0(data_dir,'/*/*-sample-IDs'))),
#                              pattern='cellrangerID_to_sangerID.txt', full.names=T)
# for(f in sampleID_files){
#   tmp_df <- read.table(f,header=T,sep='\t')

#   samplesheet <- rbind(samplesheet, 
#                        tmp_df[,c("sample_id","description",
#                                  "sanger_name","sample_name")])
  
# }
#####

#### Read .h5 files, save seurat objects ####
# Loop over the cellranger .h5 files. If it includes CITE-data, use dsb to correct 
# for background and isotype control levels
# for(f in c(filtered_feature_files,sample_feature_files)){

  # Exctract sample ID from file path
  # depending on if cellranger is run with cellranger count or cellranger multi
  # the directory structure is different, so extracting sample ID from directory
  # depends on which one is being used. This is tested with Cellranger version 6
  # and 7, might not work with other versions
  # cellranger_dir <- dirname(f)

  # if(basename(f) == 'filtered_feature_bc_matrix.h5'){
  #   sample_id <- basename(cellranger_dir)
  # }else if(basename(f) == 'sample_feature_bc_matrix.h5' || basename(f)=='sample_filtered_feature_bc_matrix.h5'){
  #   sample_id <- basename(dirname(cellranger_dir))
  # }else{
  #   stop(paste0('basename(f) should be sample_feature_bc_matrix.h5 or filtered_feature_bc_matrix.h5, but was ',basename(f)))
  # }
  
  # sample_name <- samplesheet[samplesheet$sample_id==sample_id,]$sample_name
  # print(sample_name)
  
  # name of intermediary output file
  cite_file <- paste0(tmp_rds_dir,'/1.CITE/',
                      sample_name,'.withADT.RDS')
  print(paste0(".",filtered_cellranger))
  cells = Seurat::Read10X(filtered_cellranger)

  # If there is only one assay length(cells) will be large (i.e. number of cells)
  # If there are multiple assays this number will be small
  # e.g. if GEX and CITE assay, length(cells)==2. Put it to 
  # 10 in case future studies have more assays than just GEX/CITE
  # If it has multiple assay, check that the name of one of the
  # assays is 'Antibody Capture' (otherwise, could be GEX + some other assay that's not CITE)
    # Use QC'ed cells from single cell rnaseq pipeline (from the scpred seurat file)
  # qced_cells <- readRDS(scpred_file_with_qc)
  qced_cells <- LoadH5Seurat(paste('tmp',"h5seurat",sep='.'),assays = "RNA")
  qced_barcodes <- gsub('_.*','',colnames(qced_cells))
  qced_barcodes <- gsub('-cellranger.*','',colnames(qced_cells))
  qced_cells = RenameCells(qced_cells, new.names = qced_barcodes)

  genes <-  read.table(paste0(cellranger_rawfile_path,"/features.tsv.gz"), sep = "\t", col.names = c("ENSG_ID","Gene_ID", "FeatureType"))
  qced_cells[["RNA"]] <- AddMetaData(qced_cells[["RNA"]], rownames(qced_cells[["RNA"]]), col.name = "ENSG_ID")
  qced_cells[["RNA"]] <- AddMetaData(qced_cells[["RNA"]], genes[match(rownames(qced_cells[["RNA"]]),genes$ENSG_ID),]$Gene_ID, col.name = "Gene_ID")
  cts = qced_cells@assays$RNA@counts
  rownames(cts) = genes[match(rownames(cts),genes$ENSG_ID),]$Gene_ID
  qced_cells <- CreateSeuratObject(cts, meta.data = qced_cells[[]])



  saveRDS(qced_cells, file=cite_file)

  if(length(cells)>10 | !'Antibody Capture' %in% names(cells)){
    print(paste0(sample_name,' does not have CITE assay, skip'))
    quit() 
  }
  print(paste0(sample_name,' has CITE assay, start processing'))
  
  # For dsb, we need all droplets, not only droplets that cellranger selected
  # as cells. Get the filepath of .h5 file that contains all droplets, and 
  # read into raw
  # cellranger count and cellranger multi give different directory structures, so need to
  # read raw cells in differently. cellranger multi has sample_feature_bc_matrix.h5 
  # with cells, cellranger count has filtered_feature_bc_matrix.h5 with cells
  # if(basename(f) == 'filtered_feature_bc_matrix.h5'){
  #   raw = Seurat::Read10X_h5(paste0(cellranger_dir,'/raw_feature_bc_matrix.h5'))  
  # }else if(basename(f)=='sample_feature_bc_matrix.h5' || basename(f) == 'sample_filtered_feature_bc_matrix.h5'){
    # with cellranger multi the raw file is 3 directories higher than filtered file, then in multi/count
  rna = raw <- Read10X(cellranger_rawfile_path)
  # }else{
  #   stop(paste('basename(f) should be sample_feature_bc_matrix.h5 or filtered_feature_bc_matrix.h5, but was:',basename(f)))
  # }
  
  # define cell-containing barcodes and separate cells and empty drops
  stained_cells = colnames(cells$`Gene Expression`)
  background = setdiff(colnames(raw$`Gene Expression`),stained_cells)
  
  # split the data into separate matrices for RNA and ADT
  prot = raw$`Antibody Capture`
  rna = raw$`Gene Expression`
  
  # Calculate mt % per cell, used for plotting only
  mtgene = grep(pattern = "^MT-", rownames(rna), value = TRUE) # used below
  
  # create metadata of droplet QC stats used in standard scRNAseq processing
  md = data.frame(
    rna.size = log10(Matrix::colSums(rna)), 
    prot.size = log10(Matrix::colSums(prot)), 
    n.gene = Matrix::colSums(rna > 0), 
    mt.prop = Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
  )
  
  # add indicator for barcodes Cell Ranger called as cells
  md$drop.class = ifelse(rownames(md) %in% stained_cells, 'cell', 'background')
  
  # remove barcodes with no evidence of capture in the experiment
  md = md[md$rna.size > 0 & md$prot.size > 0, ]
  
  # Set minimum and maximum protein levels and rnaseq level that will be used as
  # background. Any droplet with
  #   min_prot_size < droplet < max_prot_size & droplet < rna_size
  # will be included as background. Check below plot that this makes sense
  min_prot_size = 1.5
  max_prot_size = 3
  min_rna_size = 2.5
  ggplot(md, aes(prot.size, rna.size, colour=mt.prop))+
    geom_point(size=0.5)+
    facet_wrap(~drop.class)+
    scale_colour_viridis()+
    geom_segment(aes(x=min_prot_size,xend=max_prot_size,y=rna.size,yend=rna.size),
                 lty=2, colour='red')+
    geom_segment(aes(x=min_prot_size,xend=max_prot_size,y=0,yend=0),
                 lty=2, colour='red')+
    geom_segment(aes(x=min_prot_size,xend=min_prot_size,y=0,yend=rna.size),
                 lty=2, colour='red')+
    geom_segment(aes(x=max_prot_size,xend=max_prot_size,y=0,yend=rna.size),
                 lty=2, colour='red')+
    theme_pubr()+
    xlab('log10(total protein count)')+
    ylab('log10(total RNA count')+
    theme(legend.position='right')
  ggsave(paste0(cite_data_dir,sample_name,'.background-vs-cell.pdf'),
         width=6, height=3)

  # Set the background drops with above thresholds
  background_drops = rownames(
    md[ md$prot.size > min_prot_size & 
          md$prot.size < max_prot_size & 
          md$rna.size < min_rna_size, ]
  ) 
  background.adt.mtx = as.matrix(prot[ , background_drops])
  

  # qced_barcodes <- paste0(qced_barcodes, '_',sample_id)
  # Select the QC'ed cells from all the cellranger cells
  # colnames(prot) <- gsub('-1','',colnames(prot))
  prot_qced_cells = as.matrix(prot[ , qced_barcodes])

  # Get protein names
  prot_names <- rownames(background.adt.mtx)
  
  # define isotype controls 
  isotype.controls = prot_names[grepl('isotype',prot_names, ignore.case = T)]
  
  # normalize and denoise with dsb, use isotype controls and background droplets
  # Set return stats to get statistics for each protein
  cells.dsb.norm = DSBNormalizeProtein(
    cell_protein_matrix = prot_qced_cells, 
    empty_drop_matrix = background.adt.mtx, 
    denoise.counts = TRUE, 
    use.isotype.control = TRUE, 
    isotype.control.name.vec = isotype.controls,
    return.stats = T
  )
  
  # Save the protein stats to files for later analysis
  saveRDS(cells.dsb.norm$technical_stats, file=paste0(cite_data_dir,'/',sample_name,'.technical_stats.RDS'))
  saveRDS(cells.dsb.norm$protein_stats, file=paste0(cite_data_dir,'/',sample_name,'.protein_stats.RDS'))

  # barcoded in qced_cells object contain the sample id, add this to the dsb barcodes also. same for raw matrix
  # colnames(cells.dsb.norm$dsb_normalized_matrix) <- paste0(colnames(cells.dsb.norm$dsb_normalized_matrix), 
  #                                                          '_',sample_id)
  # colnames(prot_qced_cells) <- paste0(colnames(prot_qced_cells), 
  #                                                          '_',sample_id)
  
  # Add the two CITE assays to the QC'd cells seurat object. Note that CITE_bgRemoved is
  # added to data slot, not counts slot, because it is already normalized.
  # See https://cran.r-project.org/web/packages/dsb/vignettes/end_to_end_workflow.html
  # for more info
  qced_cells[["CITE_bgRemoved"]] = CreateAssayObject(data = cells.dsb.norm$dsb_normalized_matrix)
  qced_cells[["CITE"]] <- CreateAssayObject(counts = prot_qced_cells)
  
  current_assay <- qced_cells@active.assay

  DefaultAssay(qced_cells) <- 'CITE'
  qced_cells <- qced_cells %>% 
    NormalizeData(normalization.method = "CLR", margin = 2, verbose=F) %>%
    ScaleData(verbose=F) %>% FindVariableFeatures(verbose=F) %>% RunPCA(verbose=F)
  DefaultAssay(qced_cells) <- current_assay

  saveRDS(qced_cells, file=cite_file)
