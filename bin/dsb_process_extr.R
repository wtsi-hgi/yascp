#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)
library(Seurat)
library(tidyverse)
library(sctransform)
library(SeuratDisk)
library(ggrastr)
library(dsb)


    cellranger_filepath <- args[1] #'/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/ania/analysis_trego/work/6e/ffddc30b501cb3d43c5e76eef66715/cellbender_FPR_0.1_filtered.h5'
    # cellranger_filepath <-'cellbender_FPR_0.1_filtered.h5'
    #raw data for dbs
    cellranger_rawfile_path <- args[2] #'/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/ania/analysis_trego/work/6e/ffddc30b501cb3d43c5e76eef66715/CTRL_D1_BM__gex_data'
    # cellranger_rawfile_path <- 'CTLA4_A1_BM__gex_data'
    #AB data
    ab_filepath <- args[3] #'/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/ania/analysis_trego/work/6e/ffddc30b501cb3d43c5e76eef66715/CTRL_D1_BM__ab_data'
    # ab_filepath <- 'CTLA4_A1_BM__ab_data'
    
  #main GEX/CITEseq data
    sample <- args[4]
    # sample = 'CTLA4_A1_BM'

    ## --------------------------------------------------------------------------------------------------------------------------------------
    # We are loading the protein data, all the rna data and the cellbender background removed gex matrix
    rna = raw <- Read10X(cellranger_rawfile_path)
    cells <- Read10X_h5(cellranger_filepath, use.names = TRUE, unique.features = TRUE)
    prot = antibody <- Read10X(ab_filepath)

    # We load the 
    print("Defining bg and fg..")
    stained_cells <- colnames(cells)
    background <- setdiff(colnames(raw), colnames(cells))

    # create metadata of droplet QC stats used in standard scRNAseq processing
    # md data conmtains the mitochondrial percentages, protein counts, gene counts
    
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
    # lower_prot.size=1.5
    # upper_prot.size=2.5
    # min_rna.size=2.5

    cellmd <- md[md$drop.class == 'cell', ]
    rna.mult = (3*mad(cellmd$rna.size))
    prot.mult = (3*mad(cellmd$prot.size))
    min_rna.size = rna.lower = median(cellmd$rna.size) - rna.mult
    rna.upper = median(cellmd$rna.size) + rna.mult
    lower_prot.size = prot.lower = median(cellmd$prot.size) - prot.mult
    upper_prot.size = prot.upper = median(cellmd$prot.size) + prot.mult

    g <- md%>%
    ggplot(aes(prot.size, rna.size, colour=drop.class))+
    ggrastr::geom_point_rast(size=0.5, alpha=0.3)+
    # scale_colour_viridis()+
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
    qc_cells = rownames(
      cellmd[cellmd$prot.size > prot.lower & 
            cellmd$prot.size < prot.upper & 
            cellmd$rna.size > rna.lower & 
            cellmd$rna.size < rna.upper & 
            cellmd$mt.prop < 0.1, ]
    )

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

    stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.adt.raw))))
    stopifnot(isTRUE(all.equal(rownames(cellmd), colnames(cell.rna.raw))))

    s = Seurat::CreateSeuratObject(counts = cell.rna.raw, 
                                  meta.data = cellmd,
                                  assay = "RNA", 
                                  min.cells = 20)
    s[["CITE"]] = Seurat::CreateAssayObject(data = cells.dsb.norm$dsb_normalized_matrix,names.delim = "_")

    SaveH5Seurat(s, filename=paste0(sample,'_firstrun_dsb.h5Seurat'), overwrite = TRUE)

    # check_if_loads <- LoadH5Seurat(paste0(sample,'_firstrun_dsb.h5Seurat'))