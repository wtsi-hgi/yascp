#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
library(Seurat)
library(tidyverse)
library(sctransform)
library(SeuratDisk)
library(ggrastr)
library(dsb)
## ----sample_declaration----------------------------------------------------------------------------------------------------------------
#experiment-specific part where I get sample name and info etc.
## --------------------------------------------------------------------------------------------------------------------------------------


## --------------------------------------------------------------------------------------------------------------------------------------
#input: a clean Seurat object (to extract barcode names and to slide the dsb output in an empty slot)
# location of cellranger output: raw and filtered
# cellranger_rawfile_path raw_feature_bc_matrix.h5
# cellranger_filepath sample_filtered_feature_bc_matrix.h5

# sample.seurat <- LoadH5Seurat(file.path(output_dir, paste0(sample,"_firstrun.h5Seurat")))
# qc_cells <- colnames(sample.seurat)

#main GEX/CITEseq data
cellranger_filepath <- args[1] #'/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/ania/analysis_trego/work/6e/ffddc30b501cb3d43c5e76eef66715/cellbender_FPR_0.1_filtered.h5'
# cellranger_filepath <-'cellbender_FPR_0.1_filtered.h5'
#raw data for dbs
cellranger_rawfile_path <- args[2] #'/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/ania/analysis_trego/work/6e/ffddc30b501cb3d43c5e76eef66715/CTRL_D1_BM__gex_data'
# cellranger_rawfile_path <- 'CTRL_B1_T__gex_data'
#AB data
ab_filepath <- args[3] #'/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/ania/analysis_trego/work/6e/ffddc30b501cb3d43c5e76eef66715/CTRL_D1_BM__ab_data'
# ab_filepath <- 'CTRL_B1_T__ab_data'
# sample = strsplit(ab_filepath,'/')
# sample = sample[[1]][length(sample[[1]])]
# sample = 'CTLA4_C1_T'
# cellranger_filepath <-'cellbender_FPR_0.1_filtered.h5'
# cellranger_rawfile_path <- 'CTLA4_C1_T__gex_data'
# ab_filepath <- 'CTLA4_C1_T__ab_data'

sample <- args[4]
# sample = 'CTRL_B1_T'
## --------------------------------------------------------------------------------------------------------------------------------------

  rna = raw <- Read10X(cellranger_rawfile_path)
  cells <- Read10X_h5(cellranger_filepath, use.names = TRUE, unique.features = TRUE)
  prot = antibody <- Read10X(ab_filepath)

  
  # Convert(ab_filepath, ".h5seurat", overwrite = TRUE)
  # pbmc3k <- LoadH5Seurat("antibody-CTLA4_B1_BM.h5seurat")
  # antibody <- LoadH5Seurat("/lustre/scratch123/hgi/mdt2/teams/hgi/mo11/tmp_projects/ania/analysis_trego/work/7f/fed4885ad6ae9f4eef8142be635972/antibody-CTLA4_C1_BM.h5Seurat")
  # Convert("pbmc3k_final.h5ad", dest = "h5seurat", overwrite = TRUE)
  # antibody <- Read10X_h5('/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/ania/analysis_trego/results_rmCiteseq2/citeseq/CTLA4_C1_BM/antibody-CTLA4_C1_BM.h5ad')
  # define cell-containing barcodes and separate cells and empty drops
  print("Defining bg and fg..")
  stained_cells <- colnames(cells)
  background <- setdiff(colnames(raw), colnames(cells))

  # split the RAW data into separate matrices for RNA and ADT
  #renaming by correct_abs_names() - needs a separate file with new names (I remap antibody names to more human readable names)

  
  
  # create metadata of droplet QC stats used in standard scRNAseq processing
  #  mtgene # used below for plotting only
  
  mtgene <- grep(pattern = "^MT-", rownames(rna), value = TRUE)
  md <- data.frame(
    rna.size = log10(Matrix::colSums(rna)), 
    prot.size = log10(Matrix::colSums(prot)), 
    n.gene = Matrix::colSums(rna > 0), 
    mt.prop = Matrix::colSums(rna[mtgene, ]) / Matrix::colSums(rna)
  )
  
  # add indicator for barcodes Cell Ranger called as cells
  md$drop.class <- ifelse(rownames(md) %in% stained_cells, 'cell', 'background')
  
  # remove barcodes with no evidence of capture in the experiment (have to have prot AND rna)
  md <- md[md$rna.size > 0 & md$prot.size > 0, ]
  
  # suggested cuttofs: md$prot.size >2.5, md$rna.size < 2.5 , md$rna.size>1.8  for mt-high content droplets
  # md$prot.size >1.5, md$prot.size < 2.5 , md$rna.size> 1.8  for mt-high content droplets
  #
  #at the moment these are hard-coded; the plot produced below allows assessment whether thse cutoffs are sensible or need tweaking
    lower_prot.size=1.5
    upper_prot.size=2.5
    min_rna.size=2.5

 

    g <- md%>%
    ggplot(aes(prot.size, rna.size, colour=drop.class))+
    ggrastr::geom_point_rast(size=0.5, alpha=0.3)+
 #   scale_colour_viridis()+
    geom_segment(aes(x=lower_prot.size,xend=upper_prot.size, y=min_rna.size, yend=min_rna.size),
                 lty=2, colour='red')+
    geom_segment(aes(x=lower_prot.size,xend=upper_prot.size,y=0,yend=0),
                 lty=2, colour='red')+
    geom_segment(aes(x=lower_prot.size, xend=lower_prot.size,y=0,yend=min_rna.size),
                 lty=2, colour='red')+
    geom_segment(aes(x=upper_prot.size,xend=upper_prot.size,y=0,yend=min_rna.size),
                 lty=2, colour='red')+
    xlab('log10(total protein count)')+
    ylab('log10(total RNA count')+
    theme(legend.position='right')+theme_bw()+guides( colour = guide_legend("Cells by cellbender"))

    
#outputs_plot_with_cutoffs      
  ggsave(plot = g,paste0(sample,'.background-vs-cell_dsb.pdf'),
         width=6, height=3)

  background_drops <- rownames(
    md[ md$prot.size > lower_prot.size & 
          md$prot.size < upper_prot.size & 
          md$rna.size < min_rna.size, ]
  ) 
  
  ###Here background is defined
  background.adt.mtx <- as.matrix(prot[ , background_drops])
  
  
  #####
  # qc_cells <- colnames(cells)
  cellmd <- md[md$drop.class == 'cell', ]
  
  rna.mult = (3*mad(cellmd$rna.size))
  prot.mult = (3*mad(cellmd$prot.size))
  rna.lower = median(cellmd$rna.size) - rna.mult
  rna.upper = median(cellmd$rna.size) + rna.mult
  prot.lower = median(cellmd$prot.size) - prot.mult
  prot.upper = median(cellmd$prot.size) + prot.mult

  qc_cells = rownames(
    cellmd[cellmd$prot.size > prot.lower & 
          cellmd$prot.size < prot.upper & 
          cellmd$rna.size > rna.lower & 
          cellmd$rna.size < rna.upper & 
          cellmd$mt.prop < 0.14, ]
    )


  # cell.adt.raw <- as.matrix(prot)
  
  cell.adt.raw <- as.matrix(prot[ , qc_cells])
  cell.rna.raw = rna[ ,qc_cells]
  cellmd = cellmd[qc_cells, ]
  
  
  isotype.controls <- grep("sotype", rownames(cell.adt.raw), val=T)
  print("Normalising and denoising..")
  # normalize and denoise with dsb with default approach
  
  cells.dsb.norm <- DSBNormalizeProtein(
    cell_protein_matrix = cell.adt.raw, 
    empty_drop_matrix = background.adt.mtx, 
    denoise.counts = TRUE, 
    use.isotype.control = TRUE, 
    isotype.control.name.vec = isotype.controls,
    return.stats = T
  )
  
  #outputs_stats  
  saveRDS(cells.dsb.norm$technical_stats, file=paste0(sample,'.dsb_technical_stats.RDS'))
  saveRDS(cells.dsb.norm$protein_stats, file=paste0(sample,'.dsb_protein_stats.RDS'))

  #outputs_seurat: insert into a new slot, save 

stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.adt.raw))))
stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.rna.raw))))



# create Seurat object note: min.cells is a gene filter, not a cell filter
# s = CreateSeuratObject(counts = cell.rna.raw, 
#                                meta.data = cellmd,
#                                assay = "RNA", 
#                                )
s = Seurat::CreateSeuratObject(counts = cell.rna.raw, 
                               meta.data = cellmd,
                               assay = "RNA", 
                               min.cells = 20)
s[["CITE"]] = Seurat::CreateAssayObject(data = cells.dsb.norm$dsb_normalized_matrix,names.delim = "_")

# prots = rownames(s@assays$CITE@data)[1:28]
# s = Seurat::FindNeighbors(object = s, dims = NULL,assay = 'CITE', 
#                           features = prots, k.param = 30, 
#                           verbose = FALSE)
# c2[["ADT_dbs"]] = CreateAssayObject(data = cells.dsb.norm$dsb_normalized_matrix,names.delim = "_",)
SaveH5Seurat(s, filename=paste0(sample,'_firstrun_dsb.h5Seurat'), overwrite = TRUE)

# check_if_loads <- LoadH5Seurat(paste0(sample,'_firstrun_dsb.h5Seurat'))