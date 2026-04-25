#!/usr/bin/env Rscript
# 1. integrate all donors by SCT to get one integrated SCT object
# 2a. integrate all donors by CITE to get one integrated CITE object
# 2b. integrate all donors by CITE_bgRemoved to get one integrated CITE_bgRemoved object
# 3a. integrate SCT+CITE
# 3b. integrate SCT+CITE_bgRemoved

#### libraries ####
library("Matrix")
library(Seurat)
library(ggplot2)
library("RColorBrewer")
library(future)
options(future.globals.maxSize= 1020971520000)
if (future::supportsMulticore()) {
  future::plan(future::multicore)
} else {
  future::plan(future::multisession)
}
future.seed=TRUE
args = commandArgs(trailingOnly=TRUE)

# args=vector(mode='list', length=6); args[[1]]='pct_counts_gene_group__mito_transcript'; args[[2]]=5; args[[3]]=30; args[[4]]=12;args[[5]]= 10;args[[6]]= 9
if (args[1]=='NONE') {
  vars_to_regress = c()
  reg_name = 'regress__NONE'
}else{
  vars_to_regress = strsplit(as.character(args[1]), split=";")
  vars_to_regress = unlist(vars_to_regress)
  reg_name = paste0('regress__',args[1])
}
print(vars_to_regress)
# integrate sct
k.anchor=strtoi(args[2])
print(k.anchor)
dims=strtoi(args[3])

# Number of dimensions to use for SCT integrated data (check elbow plot to see if it's correct)
ndim_sct = strtoi(args[4])

# Number of dimensions to use for cite_bgRemoved integrated data (check elbow plot to see if it's correct)
ndim_citeBgRemoved = strtoi(args[5])

# Number of dimensions to use for cite integrated data (check elbow plot to see if it's correct)
ndim_cite_integrated = strtoi(args[6])
 

#####
# .libPaths('/software/hgi/containers/yascp/jaguar_rlibs')
# install.packages("Matrix")
#### set up directories, colors paths ####
data_dir <- getwd()
outdir <- paste0('.')
# dir.create(outdir,showWarnings = F)
figdir <- paste0(outdir,'/figures','__',reg_name,'/')
dir.create(figdir,showWarnings = F)
# tmp_rds_dir <- paste0(outdir,'/tmp_rds_files/')
# dir.create(tmp_rds_dir,showWarnings = F)
tmp_rds_file <- paste0(reg_name,'__','all_samples_integrated.RDS')

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

# Select the sc eqtlgen rds files
wg2_dir <- getwd()


cite_files <- list.files(pattern='.RDS',
                         path=Sys.glob('.'),recursive=T, full.names=T)
#####

#### Get samplesheets for all batches and extract sample IDs ####
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

##### Get the donor info per cell from vireo and from matchings we did #####
# read in donor info so that integration can be done per donor
print('---- Combine files together -----')
# Get the location of the donor files
donor_files <- list.files(pattern='donor_ids.tsv',
                          path=Sys.glob(paste0(data_dir,'/vireo_*')),
                          full.names=T)

matched_donor_files <- 'matched_donors.txt'
matched_donors <- read.table(matched_donor_files,header=T,sep='\t')

donor_cells <- data.frame()
count = 1
drop_doublets = TRUE
for(f in donor_files){
  sample_id <- basename(dirname(f))

  sample_name <-gsub('vireo_','',sample_id)
  df <- read.table(f, header=T, sep='\t')

  df$sample_name <- sample_name
  df$sample <- sample_name
  
  # This is still some old code from when there were donors in multiple batches but when we didn't
  # have genotypes, at some point his will be updated. Get the matched genotypes
  
  df$doublet_logLikRatio <- NULL
  if (drop_doublets){
    print('--- Dropping Doublets ---')
    df <- df[!grepl('doublet|unassigned', df$donor_id),]
  }else{
    print('--- Not Dropping Doublets ---')
    df <- df[!grepl('unassigned', df$donor_id),]
    matched_donors[nrow(matched_donors) + 1,] = c(paste0("doublet",count),sample_name,"doublet")
  }
  
  df$matched.donor <- matched_donors[match(paste0(df$sample_name, 
                                                  df$donor_id), 
                                           paste0(matched_donors$sample, 
                                                  matched_donors$old_donor)),]$new_donor
  donor_cells <- rbind(donor_cells, df)
  count = count+1
}
#####

##### SCTransforming and correcting RNA-seq data, calculate cellcycle, remove MT/ribosomal genes #####
print('---- Normalise and prepeare -----')
sobj_list <- list()
sample_names <- c()
i <- 1
for(f in cite_files){
  sample_id <- gsub('.withADT|.RDS','',basename(f))
  sample_name <- sample_id

  #### for testing only pool 1-5 ####
  #if(!grepl('01|02|03|04|05', sample_name)){
  #  next
  #}
  #####
  print(sample_name)

  # for samples that have CITE-seq data, use sobj with CITE data included (from 1.add_ADT.R)
  # else, use the scpred step output file.
  # if(sum(grepl(sample_name,cite_files))== 0){
  # print(paste(sample_name,' does not have cite data, use scpred file'))
  sobj <- readRDS(f)  #// Here we actually dont have the doublets anymore since they were not merged back in the files.
  # t2 = readRDS('/lustre/scratch123/hgi/teams/hgi/mo11/tmp_projects/jaguar_yascp/nieks_pipeline/sle-niek/2022-03-07-sc-eQTLgen-pipeline/2022-03-17-SLE-datasets/2022-07-13-SLE/no-cellbender/WG1/2023-11-23-batch20/2.Demultiplexing-and-Doublet-Removal/out/QC_figures/seurat_object_all_pools_all_barcodes_final_assignments.rds')
  # }else if(sum(grepl(sample_name,cite_files))==1){
  #   print(paste(sample_name,'has cite data, use cite file'))
  #   cfile <- cite_files[grepl(sample_name,cite_files)]
  #   sobj <- readRDS(cfile)  
  # }else{
  #   # Just a safety check, should never get here
  #   stop('Grepping smample name should give one file, but more found. Check code')
  # }
  
  # remove mito and ribo genes
  # Selecting @assays$RNA because if using subset() or 
  # sobj[!grepl("^MT-", rownames(sobj)), ], the CITE-seq assays get removed from object
  sobj@assays$RNA <- sobj[!grepl("^MT-", rownames(sobj)), ]@assays$RNA
  sobj@assays$RNA <- sobj[!grepl("^RP[SL]", rownames(sobj)), ]@assays$RNA
  sobj$Pool <- sample_name
  

  # Add the vireo donor assignment to the seurat object
  this_donor_cells <- donor_cells[donor_cells$sample==sample_id,]

  if (dim(this_donor_cells)[1]==0){
    # Here we skip non deconvoluted samples
    next
  }
  sobj@meta.data$donor.vireo <- this_donor_cells[match(paste0(rownames(sobj@meta.data),'_',sample_name),
                                                       paste0(this_donor_cells$cell,
                                                              '_',this_donor_cells$sample)),]$matched.donor
    # sobj@meta.data$donor.vireo <- this_donor_cells[match(rownames(sobj@meta.data),
    #                                                    this_donor_cells$cell),]$matched.donor
  # remove NA values as these are either doublet or unassigned
  sobj <- sobj[,rownames(sobj@meta.data[!is.na(sobj$donor.vireo),])]

  # make one seurat object per donor for this pool  
  sobj_per_donor <- SplitObject(sobj, split.by = "donor.vireo")
  for(donor in names(sobj_per_donor)){
    print(donor)
  
    # Regress only percent mito,
    # because sctransform should already control for read count and n genes
    # Using glmGamPoi because faster and more robust (see also https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02584-9)
    # return all genes so that more can be used for calculating cell cycle score
    sobj_per_donor[[donor]] <-SCTransform(sobj_per_donor[[donor]], 
                        vars.to.regress=vars_to_regress, 
                        method = "glmGamPoi",
                        verbose = F, return.only.var.genes = F)

    # Calculate cell cycle scores
    sobj_per_donor[[donor]] <- CellCycleScoring(
      NormalizeData(sobj_per_donor[[donor]]),
      s.features = cc.genes$s.genes,
      g2m.features = cc.genes$g2m.genes,
      assay = 'SCT',
      set.ident = TRUE
    )
    
    # Calculate and plot PCA on cell cycle genes only first
    sobj_per_donor[[donor]] <-  RunPCA(sobj_per_donor[[donor]], features = c(cc.genes$s.genes, cc.genes$g2m.genes),verbose=F)
    p1 <- DimPlot(sobj_per_donor[[donor]])+ggtitle('Before cell cycle correction (cell cycle genes only)')
    # Then also on all variable genes
    sobj_per_donor[[donor]] <-  RunPCA(sobj_per_donor[[donor]], dims=1:10,verbose=F)
    p2 <- DimPlot(sobj_per_donor[[donor]])+ggtitle('Before cell cycle correction (all genes)')
    
    # Would regress out cell cycle here, but most of the time the cell cycle genes
    # are mostly not expressed in dataset (at least in batch 1), which makes me
    # not trust the scores. So, only regress mitochondria
    sobj_per_donor[[donor]] <- SCTransform(sobj_per_donor[[donor]], 
                        vars.to.regress=vars_to_regress,
                                          #"S.Score", "G2M.Score"), # see comments above, don't want to regress cell cycle score
                        method = "glmGamPoi",
                        verbose = F)
    
    #### below block is plotting PCs after cell cycle regression, which showed little effect ####
    # Since cell cycle score is not regressed anymore, this plotting no longer needed
    # Calculate and plot PCA on cell cycle genes only first (after cell cycle regression)
    #sobj_per_donor[[donor]] <-  RunPCA(sobj_per_donor[[donor]], features = c(cc.genes$s.genes, cc.genes$g2m.genes), verbose=F)
    #p3 <- DimPlot(sobj_per_donor[[donor]])+ggtitle('After cell cycle correction (cell cycle genes only)')
    # Then also on all variable genes
    #sobj_per_donor[[donor]] <-  RunPCA(sobj_per_donor[[donor]], dims=1:10, verbose=F)
    #p4 <- DimPlot(sobj_per_donor[[donor]])+ggtitle('Before cell cycle correction (all genes)')
    
    # Put the plots together and save
    #(p1+p2)/(p3+p4)
    #dir.create(paste0(figdir,'/2.cellCycleCorrection/'),showWarnings = F)
    #ggsave(paste0(figdir,'/2.cellCycleCorrection/',donor,'.pdf'), width=6, height=6)
    #print(paste('Figure saved to',paste0(figdir,'/2.cellCycleCorrection/',donor,'.pdf')))
    #####
    
    sobj_per_donor[[donor]]@project.name <- paste0(sample_name,'-',donor)
    # Then add this pool-donor seurat object to list so that sobj_list
    # will contain seurat objects for all donors of all pools
    sobj_list[[i]] <- sobj_per_donor[[donor]]
    sample_names <- c(sample_names,  paste0(sample_name,'-',donor))
    i <- i+1
  }
}


print("---- Integrate ------")
names(sobj_list) <- sample_names
#####
Csparse_validate = "CsparseMatrix_validate"
#### Functions for integrating sct and cite-seq data ####

integrate_sct <- function(slist, reference_samples, k.anchor=5, dims=30){
  # Integrate a list of seurat objects together using the SCT assay using rpca
  # slist: list of seurat objects
  # reference_samples: list of names to use as reference (reference parameter in FindIntegrationAnchors)
  # k.anchor: Number of anchors to use in FindIntegrationAnchors
  # dims: Number of dimeonsions to use in FindIntegrationAnchors
  # See https://satijalab.org/seurat/articles/integration_introduction.html and
  #     https://satijalab.org/seurat/articles/integration_rpca.html for more info
  
  # In case SCT is not default assay, change it to SCT
  for(i in 1:length(slist)){
    DefaultAssay(slist[[i]]) <- 'SCT'
  }

  print('Using reference samples:')
  print(reference_samples)
  
  print('Select integration features')
  features <- SelectIntegrationFeatures(object.list = slist, nfeatures = 3000)
  
  print('Prep SCT integration')
  slist <- PrepSCTIntegration(object.list = slist, anchor.features = features)
  
  print('Run PCA')
  slist <- lapply(X = slist, FUN = RunPCA, features = features,verbose=F)
  
  # reference parameter of FindIntegrationAnchors expects indexes, not names, so
  # look up the index of the samples given to reference_samples
  reference_sample_index <- which(names(slist) %in% reference_samples)
  print('Find anchors')
  # Use rpca because it is faster for large dataset, see https://satijalab.org/seurat/articles/integration_rpca.html
  Csparse_validate = "CsparseMatrix_validate"
  anchors <- FindIntegrationAnchors(object.list = slist, normalization.method = "SCT",
                                    anchor.features = features, dims = 1:dims, 
                                    k.anchor = k.anchor, reduction = 'rpca',
                                    reference=reference_sample_index)
  print('Integrate data')
  integrated.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:dims)
  return(integrated.sct)
}

integrate_cite <- function(slist,reference_samples, normalize, assay, k.anchor=5, dims=20){


  # slist =sobj_list;random_donors_for_reference=reference_samples;normalize=F; assay='CITE_bgRemoved';k.anchor=5; dims=20
                                                     
  # Integrate a list of seurat objects together using assay (e.g. CITE or CITE_bgRemoved)
  # slist: list of seurat objects
  # reference_samples: list of names to use as reference (reference parameter in FindIntegrationAnchors)
  # normalize: Whether to normalize CITE data or not
  # assay: Which assay to use
  # k.anchor: Number of anchors to use in FindIntegrationAnchors
  # dims: Number of dimeonsions to use in FindIntegrationAnchors. Lower than default because less proteins
  
  print('Using reference samples:')
  print(reference_samples)
  
  for(i in 1:length(slist)){

    DefaultAssay(slist[[i]]) <- assay
  }
  
  # Normalize date if paramter is set T

  slist_test <- lapply(X = slist, FUN = function(x) {
    if(normalize==T){
      x <- NormalizeData(x, verbose = FALSE, normalization.method='CLR', margin=2) 
    }
  })
  
  
  print('Select integration features')
  seurat.features.adt <- rownames(slist[[1]])
  
  print(paste('# of selected features before isotype removal:', length(seurat.features.adt)))
  # remove isotype controls from features as these should not be used for integration
  seurat.features.adt <- seurat.features.adt[!grepl('isotype',seurat.features.adt,
                                                    ignore.case = T)]
  print(paste('# of selected features after isotype removal:', length(seurat.features.adt)))
  print('Scale data and run PCA')
  slist <- lapply(X = slist, FUN = function(x) {
    x <- ScaleData(x, features = seurat.features.adt, verbose = FALSE)
    x <- RunPCA(x, features = seurat.features.adt, verbose = FALSE)
  })
  
  print('Find anchors')
  anchors.adt <- FindIntegrationAnchors(object.list = slist, 
                                        reduction = "rpca", 
                                        dims = 1:dims)
  
  print('Integrate data')
  integrated.adt <- IntegrateData(anchorset = anchors.adt, dims = 1:dims)
  
  # NOTE: With rpca integration for ADT we do scale after integration
  # unlike SCT integration for expression above
  print('Scale data and run PCA and umap')
  integrated.adt <- ScaleData(integrated.adt)
  return(integrated.adt)
}
#####

pools_to_integrate <- unique(gsub('-donor[0-9]+','',names(sobj_list)))

# from each pool randomly select one donor to use as reference with one donor from all other pools
random_donors_for_reference <- c()
pools_to_use <- c()

# Set a seed so that when rerunning the code the same samples get sampled
set.seed(94235)

# For each pool, randomly sample 1 donor
for(pool in pools_to_integrate){
  random_donor <- sample(names(sobj_list)[grepl(pool, names(sobj_list))],1)
  random_donors_for_reference <- c(random_donors_for_reference, random_donor)
}

# integrate sct
# slist = sobj_list
# reference_samples = random_donors_for_reference
print('---- Integrating ---')

random_donor_integration_sct <- integrate_sct(sobj_list, random_donors_for_reference)


# integrate CITE
random_donor_integration_cite <- integrate_cite(sobj_list, 
                                              random_donors_for_reference,
                                              normalize=F,
                                              assay='CITE')


# integrate cite_bgRemoved
random_donor_integration_citeBgRemoved <- integrate_cite(sobj_list, 
                                                         random_donors_for_reference,
                                                         normalize=F,
                                                         assay='CITE_bgRemoved')
#####
print("---- Integration done, plot now -----")
#### PCA and UMAP ####
dir.create(paste0(figdir,'/2.elbow_plots/'),showWarnings = F)

random_donor_integration_sct <- RunPCA(random_donor_integration_sct, verbose = FALSE)

# Number of dimensions to use for SCT integrated data (check elbow plot to see if it's correct)
ElbowPlot(random_donor_integration_sct, ndims=50)+geom_vline(xintercept = ndim_sct+0.5,
                                                             lty=2,colour='red')
ggsave(paste0(figdir,'/2.elbow_plots/SCT.pdf'), width=5, height=5)
random_donor_integration_sct <- RunUMAP(random_donor_integration_sct, reduction = "pca", dims = 1:ndim_sct)

# Number of dimensions to use for cite_bgRemoved integrated data (check elbow plot to see if it's correct)
random_donor_integration_citeBgRemoved <- RunPCA(random_donor_integration_citeBgRemoved, verbose = FALSE)
ElbowPlot(random_donor_integration_citeBgRemoved, ndims=50)+geom_vline(xintercept = ndim_citeBgRemoved+0.5,
                                                                       lty=2,colour='red')
ggsave(paste0(figdir,'/2.elbow_plots/CITE_bgRemoved.pdf'), width=14, height=4)
random_donor_integration_citeBgRemoved <- RunUMAP(random_donor_integration_citeBgRemoved, reduction = "pca", 
                                                  dims = 1:ndim_citeBgRemoved)

# Number of dimensions to use for cite integrated data (check elbow plot to see if it's correct)
random_donor_integration_cite <- RunPCA(random_donor_integration_cite, verbose = FALSE)
ElbowPlot(random_donor_integration_cite, ndims=50)+geom_vline(xintercept = ndim_cite_integrated+0.5,
                                                              lty=2,colour='red')
random_donor_integration_cite <- RunUMAP(random_donor_integration_cite, reduction = "pca", 
                                         dims = 1:ndim_cite_integrated)


#####

##### plot each umap #####
n_pools <- length(unique(random_donor_integration_sct$Pool))
n_celltype <- length(unique(random_donor_integration_sct$Azimuth.predicted.celltype.l2))

p1 <- DimPlot(random_donor_integration_sct, group.by='Pool',  cols=myPalette(n_pools), shuffle=T)
p2 <- DimPlot(random_donor_integration_sct, group.by='Azimuth.predicted.celltype.l2',  cols=myPalette(n_celltype))
p1+p2       
dir.create(paste0(figdir,'/2.initial_umaps/'),showWarnings = F)
ggsave(paste0(figdir,'/2.initial_umaps/all_samples_sct.pdf'), width=14, height=4)

p1 <- DimPlot(random_donor_integration_citeBgRemoved, group.by='Pool',  cols=myPalette(n_pools), shuffle=T)
p2 <- DimPlot(random_donor_integration_citeBgRemoved, group.by='Azimuth.predicted.celltype.l2',  cols=myPalette(n_celltype))
p1+p2       
dir.create(paste0(figdir,'/2.initial_umaps/'),showWarnings = F)
ggsave(paste0(figdir,'/2.initial_umaps/all_samples_citeBgRemoved.pdf'), width=14, height=4)

p1 <- DimPlot(random_donor_integration_cite, group.by='Pool',  cols=myPalette(n_pools), shuffle=T)
p2 <- DimPlot(random_donor_integration_cite, group.by='Azimuth.predicted.celltype.l2',  cols=myPalette(n_celltype))
p1+p2       
dir.create(paste0(figdir,'/2.initial_umaps/'),showWarnings = F)
ggsave(paste0(figdir,'/2.initial_umaps/all_samples_citeBgRemoved.pdf'), width=14, height=4)
#####


##### add the assays together #####

random_donor_integration_sct[["CITE.integrated"]] <- random_donor_integration_cite[["integrated"]]
random_donor_integration_sct[["CITE_bgRemoved.integrated"]] <- random_donor_integration_citeBgRemoved[["integrated"]]
random_donor_integration_sct[["pca.cite.integrated"]] <- random_donor_integration_cite[["pca"]]
random_donor_integration_sct[["pca.cite_bgRemoved.integrated"]] <- random_donor_integration_citeBgRemoved[["pca"]]
random_donor_integration_sct[["umap.cite.integrated"]] <- random_donor_integration_cite[["umap"]]
random_donor_integration_sct[["umap.cite_bgRemoved.integrated"]] <- random_donor_integration_citeBgRemoved[["umap"]]

  # Rename SCT assays
random_donor_integration_sct <- RenameAssays(object = random_donor_integration_sct, integrated = 'SCT.integrated')
names(random_donor_integration_sct@reductions)[[which(names(random_donor_integration_sct@reductions)=='umap')]] <- 'umap.SCT.integrated'
names(random_donor_integration_sct@reductions)[[which(names(random_donor_integration_sct@reductions)=='pca')]] <- 'pca.SCT.integrated'

# remove predicted_ADT because we have our own CITE measurements
random_donor_integration_sct@assays$predicted_ADT <- NULL
#####

print('---- Done, lets save the Object ----')
saveRDS(random_donor_integration_sct, file=tmp_rds_file)

