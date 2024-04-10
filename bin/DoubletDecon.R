#!/usr/bin/env Rscript

# .libPaths("/usr/local/lib/R/site-library")
suppressMessages(suppressWarnings(library(argparse)))


# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-o", "--out", required = TRUE, help="The output directory where results will be saved")
parser$add_argument("-s", "--seurat_object", required = TRUE, type = "character", help = "A QC, normalized seurat object with classifications/clusters as Idents() saved as an rds object.")
parser$add_argument("-g", "--num_genes", required = FALSE, type = "integer", default=50, help = "Number  of genes to use in \'Improved_Seurat_Pre_Process\' function.")
parser$add_argument("-r", "--rhop", required = FALSE, type="double", default=0.9, help="rhop to use in DoubletDecon - the number of SD from the mean to identify upper limit to blacklist")
parser$add_argument("-p", "--species", required = FALSE, type = "character", default="hsa", help = "The species of your sample. Can be scientific species name, KEGG ID, three letter species abbreviation, or NCBI ID.")
parser$add_argument("-n", "--nCores", required = FALSE, type = "double", default=2, help = "The number of unique cores you would like to use to run DoubletDecon. By default, uses one less than available detected.")
parser$add_argument("-c", "--removeCC", required = FALSE, type = "logical", default=FALSE, help = "Whether to remove clusters enriched in cell cycle genes.")
parser$add_argument("-m", "--pmf", required = FALSE, type = "logical", default=TRUE, help = "Whether to use unique gene expression in doublet determination.")
parser$add_argument("-f", "--heatmap", required = FALSE, type = "logical", default=FALSE, help = "Whether to generate heatmaps.")
parser$add_argument("-t", "--centroids", required = FALSE, type = "logical", default=FALSE, help = "Whether to use centroids instead of medoids for doublet detecting.")
parser$add_argument("-d", "--num_doubs", required = FALSE, type = "integer", default=100, help = "The number of doublets to simulate for each cluster pair.")
parser$add_argument("-5", "--only50", required = FALSE, type = "logical", default=FALSE, help = "Whether to only compute doublets as 50:50 ratio. Default is to use other ratios as well.")
parser$add_argument("-u", "--min_uniq", required = FALSE, type = "integer", default=4, help = "Minimum number of unique genes to rescue a cluster identified as doublets.")
                                        
# get command line options, if help option encountered print help and exit,
# otherwise if options not found on command line then set defaults, 
args <- parser$parse_args()

suppressMessages(suppressWarnings(library(DoubletDecon)))
suppressMessages(suppressWarnings(library(tidyverse)))
suppressMessages(suppressWarnings(library(Seurat)))
suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(data.table)))
library(viridis)
library(SeuratDisk)
library(future)
# options(future.globals.maxSize= 1020971520000)
# options(na.action="na.exclude")
# options(na.action="na.omit")
# if (future::supportsMulticore()) {
#   future::plan(future::multicore)
# } else {
#   future::plan(future::multisession)
# }

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
seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
print('PCA performed')
seurat <- FindNeighbors(seurat, dims = 1:10)
print('Neighbors found')
seurat <- FindClusters(seurat, resolution = 0.5)
print(seurat[["pca"]], dims = 1:5, nfeatures = 5)

# seurat <- Read10X('TMP_DIR')
# seurat_object = CreateSeuratObject(counts = seurat)
## Preprocess ##
processed <- Improved_Seurat_Pre_Process(seurat, num_genes=args$num_genes, write_files=FALSE)

DeconRNASeq = function(datasets, signatures, proportions=NULL, checksig=FALSE, known.prop = FALSE, use.scale = TRUE, fig = TRUE){

  if (is.null(datasets)) 
      stop(" Missing the mixture dataset, please provide a tab-delimited text file for mixture samples.")
  if (is.null(signatures)) 
      stop(" Missing the signature dataset, please provide a tab-delimited text file for pure tissue/cell types.")
  if (is.null(proportions)&&known.prop) 
      stop(" Missing the known proprotions, please provide a tab-delimited text file containing known fractions for pure tissue/cell types.")


  ## load data 

  #x.signature <- read.table(signatures, header=T, row.names=1, sep="\t")
  #x.data <- read.table(datasets, sep="\t", header=T, row.names=1)
  x.signature <- signatures
  x.data <- datasets

  #error checks
  if (is.data.frame(x.signature)==F) 
      stop("signature datasets must be a dataframe")
  if (sum(is.na(x.signature))>0)
      stop("signature data cannot have NAs. please exclude or impute missing values.")
  if (is.data.frame(x.data)==F)
      stop("mixture datasets must be a dataframe")
  if (sum(is.na(x.data))>0)
      stop("mixture data cannot have NAs. please exclude or impute missing values.")


  numofg <- nrow(x.signature)
  Numofx <- ncol(x.signature)

  #Underdetermined systems check
  if (numofg < Numofx)
      stop("The number of genes is less than the number of cell types, which means less independent equations than unknowns.")


  #PCA analysis to estimate the number of pure cell types in the mixture
  x.data.temp <- prep(x.data, scale = "none", center = TRUE)
  x.data.pca <- pca(x.data.temp, method = "svd", center = FALSE, nPcs = Numofx)
  #summary of PCA
  out.pca <- summary(x.data.pca)
  Var <- R2cum(x.data.pca)
  numofmix <- order(Var>0.99,decreasing=T)[1]
  
  if (fig) {
     if (Numofx != numofmix) {
        cat("\n Attention: the number of pure cell types =", Numofx, " defined in the signature matrix;\n" )
        cat("\n PCA results indicate that the number of cell types in the mixtures =", numofmix, "\n" )
     }
  }
  
  if (checksig && numofg >=40){
      #generate the steps to check the condition number of the signature matrix 
      step <- seq(20,numofg, by=20) #every 20 genes
      #step-wise compute the condition number for the subsets of sigature
      sig.cond <- sapply(step, function(x) kappa(scale(x.signature[1:x,]))) 
      condplot(step, sig.cond)
  }
  #x.data.pca <- prcomp(x.data, center=T, scale.=T, retx=T)
  

  ## reorder the rows to match the signature file
  common.signature <- rownames(x.signature) %in% rownames(x.data)
  common.data <- rownames(x.data) %in% rownames(x.signature)
  x.data <- x.data[common.data,]
  x.signature <- x.signature[common.signature,]
  
  x.subdata <- x.data[rownames(x.signature),]

  ## number of cell/tissue types
  Numofx <- ncol(x.signature)

  ## quadratic programming preparation
  if (use.scale) {
     AA <- scale(x.signature)
  }
  else {
     AA <- x.signature
  }

  EE <- rep(1, Numofx)
  FF <- 1
  GG <- diag(nrow=Numofx)
  HH <- rep(0, Numofx)

  out.all <- c()

  for (i in colnames(x.subdata)) {

    BB <- x.subdata[,i]
    if (use.scale) {
      if (all(BB == 0)){
        BB <- BB
      }else{
        BB <- scale(BB)
      }
    }
  
    out <- lsei(AA, BB, EE, FF, GG, HH)
    out.all <- rbind(out.all, out$X)

  }
  
  mean.rmse <- 0
  rmse <- c()
  if (known.prop)
  {
   x.proportions <- proportions

   x.proportions <- x.proportions[colnames(x.data),]

   parray <- ggplot()
   length(parray) <- ncol(out.all)
   
   for (i in 1:ncol(out.all)){
       A.error <- rmse(x.proportions[,i], out.all[,i])
       rmse <- c(rmse, A.error)
       x <- out.all[,i]
       y <- x.proportions[,i]
  
       xlabel <- paste('estimated ',colnames(x.proportions)[i])
       ylabel <- paste('actual ', colnames(x.proportions)[i])
       main <- sprintf("Scatter plot of proportions,\n RMSE = %.3f", A.error)
       
       parray[[i]] <- ggplot(data.frame(x,y), aes(x,y)) + geom_point(alpha=.3)+labs(title=main)+geom_abline(intercept=0, slope=1, colour = "red", size = 1)+ xlab(xlabel) + ylab(ylabel)
        
   }
   if (fig)
   {
   	multiplot(plotlist = parray, cols=2)
   }
   mean.rmse <- mean(rmse)
  }
  
  if (known.prop){
  	return(list(out.all = out.all, out.pca = out.pca, out.rmse = mean.rmse))
  }
  else{
  	return(list(out.all = out.all, out.pca = out.pca))
  }
  
}

# Create Line Chart for condition number

condplot = function(step,cond){
  yrange <- range(cond) 
  xrange <- c(1, length(step))

  #png('condnum.png', width=500, height=500)
  # set up the plot
  plot(xrange, yrange, type="n", xlab="Number of genes in the signature", ylab="Condition number", xaxt="n")
  colors <- rainbow(1)
  linetype <- 1
 
  # add lines
  lines(1:length(step), cond, type="b", lwd=1.5,lty=linetype, col=colors)
  
  # draw an axis on the bottom
  step.ind <- seq(1, length(step), by = 3)
  axis(1, at=step.ind,labels=step[step.ind])

  # add a title and subtitle
  title("Condition number of the signature matrix")
  #dev.off()
  
}


Is_A_Doublet<-function(data, newMedoids, groups, synthProfiles, log_file_name){
  
  #create data frame to store doublets table
  isADoublet=data.frame(matrix(ncol=4,nrow=(ncol(data)-1)))
  rownames(isADoublet)=colnames(data)[2:ncol(data)]
  rownames(newMedoids)=rownames(data)[2:nrow(data)]
  
  #run DeconRNASeq with new medoids and data
  d3 = data[2:nrow(data), 2:ncol(data)]
  results=DeconRNASeq(d3, newMedoids)
  resultsreadable=round(results$out.all*100,2)
  rownames(resultsreadable)=rownames(isADoublet) #make an easily readable results table
  
  #get average profiles for cell clusters
  averagesReal=as.data.frame(matrix(ncol=ncol(resultsreadable), nrow=length(unique(groups[,2]))))
  colnames(averagesReal)=colnames(resultsreadable)
  for(clust in 1:length(unique(groups[,2]))){
    cells=row.names(subset(groups, groups[,1]==clust))
    subsetResults=resultsreadable[row.names(resultsreadable) %in% cells,]
    averagesReal[clust,]=apply(subsetResults,2,mean)
  }
  
  #create a table with average profiles of cell clusters and synthetic combinations
  allProfiles=rbind(averagesReal, synthProfiles)
  
  #this section determines the profile with the highest correlation to the given cell and determines if it is one of the doublet profiles
  for(cell in 1:nrow(isADoublet)){
    if(ncol(resultsreadable)==2){ #If there are only 2 groups, correlation won't work, so I use minimum euclidean distance instead
      a=rbind(allProfiles, resultsreadable[cell,])
      b=as.matrix(dist(a))
      c=b[nrow(b), 1:(ncol(b)-1)]
      chosenCorrelation=c[c %in% min(c)]
      isADoublet[cell,1]=100-chosenCorrelation #100-euclidean distance
      isADoublet[cell,2]=names(chosenCorrelation)
      if(names(chosenCorrelation) %in% unique(groups[,2])){ #it is an original cluster
        isADoublet[cell,3]=FALSE
      }else{
        isADoublet[cell,3]=TRUE
      }
    }else{
      #correlations=apply(allProfiles, 1, cor, resultsreadable[cell,])
      correlations=apply(allProfiles, 1, cor, resultsreadable[cell,])
      sortCorrelations=sort(correlations, decreasing = T)[1:2]
      maxCorrelation1=which(correlations==sortCorrelations[1])
      maxCorrelation2=which(correlations==sortCorrelations[2])
      chosenCorrelation=maxCorrelation1
      isADoublet[cell,1]=correlations[chosenCorrelation]
      correlatedCluster=row.names(allProfiles)[chosenCorrelation]
      isADoublet[cell,2]=correlatedCluster
      if(chosenCorrelation>length(unique(groups[,2]))){
        isADoublet[cell,3]=TRUE
      }else{
        isADoublet[cell,3]=FALSE
      }
    }
    
  }
  
  isADoublet[,4]=groups[,2]
  
  colnames(isADoublet)=c("Distance","Cell_Types","isADoublet","Group_Cluster")
  
  cat(paste0(length(which(isADoublet$isADoublet==TRUE)),"/", nrow(isADoublet),  " possible doublets removed"), file=log_file_name, append=TRUE, sep="\n")
  
  return(list(isADoublet=isADoublet, resultsreadable=resultsreadable))
  
}


Main_Doublet_Decon<-function(rawDataFile, groupsFile, filename, location, fullDataFile=NULL, removeCC=FALSE, species="mmu", rhop=1, write=TRUE, PMF=TRUE, useFull=FALSE, heatmap=TRUE, centroids=TRUE, num_doubs=100, only50=FALSE, min_uniq=4, nCores=-1){

  #load required packages
  cat("Loading packages...", sep="\n")
  library(DeconRNASeq)
  library(gplots)
  library(plyr)
  library(MCL)
  library(clusterProfiler)
  library(mygene)
  library(tidyr)
  library(R.utils)
  library(dplyr)
  library(foreach)
  library(doParallel)
  library(stringr)

  #Set up log file
  log_file_name=paste0(location, filename,".log")
  log_con <- file(log_file_name)
  cat(paste0("filename: ",filename), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("location: ",location), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("removeCC: ",removeCC), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("species: ",species), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("rhop: ",rhop), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("write: ",write), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("PMF: ",PMF), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("useFull: ",useFull), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("heatmap: ",heatmap), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("centroids: ",centroids), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("num_doubs: ",num_doubs), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("only50: ",only50), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("min_uniq: ",min_uniq), file=log_file_name, append=TRUE, sep="\n")
  cat(paste0("nCores: ",nCores), file=log_file_name, append=TRUE, sep="\n")

  #Check variables
  if(is.character(rawDataFile)!=TRUE & is.data.frame(rawDataFile)!=TRUE){print("ERROR: rawDataFile must be a character string!")}
  if(is.character(groupsFile)!=TRUE & is.data.frame(groupsFile)!=TRUE & is.matrix(groupsFile)!=TRUE){print("ERROR: groupsFile must be a character string!")}
  if(is.character(filename)!=TRUE){print("ERROR: filename must be a character string!")}
  if(is.character(location)!=TRUE){print("ERROR: location must be a character string!")}
  if(is.character(fullDataFile)!=TRUE & is.null(fullDataFile)!=TRUE & is.data.frame(fullDataFile)!=TRUE){print("ERROR: fullDataFile must be a character string or NULL!")}
  if(is.logical(removeCC)!=TRUE){print("ERROR: removeCC must be TRUE or FALSE!")}
  if(is.character(species)!=TRUE){print("ERROR: species must be a character string!")}
  if(is.numeric(rhop)!=TRUE){print("ERROR: rhop must be numeric!")}
  if(is.logical(write)!=TRUE){print("ERROR: write must be TRUE or FALSE!")}
  if(is.logical(PMF)!=TRUE){print("ERROR: PMF must be TRUE or FALSE!")}
  if(is.logical(useFull)!=TRUE){print("ERROR: useFull must be TRUE or FALSE!")}
  if(is.logical(heatmap)!=TRUE){print("ERROR: heatmap must be TRUE or FALSE!")}
  if(is.logical(centroids)!=TRUE){print("ERROR: centroids must be TRUE or FALSE!")}
  if(is.numeric(num_doubs)!=TRUE){print("ERROR: numdoubs must be numeric!")}
  if(is.logical(only50)!=TRUE){print("ERROR: only50 must be TRUE or FALSE!")}
  if(is.numeric(min_uniq)!=TRUE){print("ERROR: min_uniq must be numeric!")}
  
  #Read in data
  cat("Reading data...", file=log_file_name, append=TRUE, sep="\n")
  cat("Reading data...", sep="\n")
  
  ICGS2_flag=F #set for checking if the input file is in ICGS2 format
  
  if(class(rawDataFile)=="character"){
    #NEW: test for ICGS2
    rawDataHeader=read.table(rawDataFile, sep="\t",header=F, row.names=1, nrows=1, stringsAsFactors = F)
    if(length(grep(":", rawDataHeader[2]))==1){
      ICGS2_flag=T
      ICGS2=ICGS2_to_ICGS1(rawDataFile, groupsFile, log_file_name)
      rawData=ICGS2$rawData
    }else{
      rawData=read.table(rawDataFile, sep="\t",header=T, row.names=1, stringsAsFactors = T)
    }
  }else{
    cat("WARNING: if using ICGS2 file input, please import 'rawDataFile' and 'groupsFile' as path/location instead of an R object." , sep="\n")
    rawData=rawDataFile
  }

  if(class(groupsFile)[1]=="character"){
    if(ICGS2_flag==T){
      groups=ICGS2$groups
    }else{
      groups=read.table(groupsFile, sep="\t",header=F, row.names=1, stringsAsFactors = T)
    }
  }else{
    groups=groupsFile
  }

  #Clean up data and groups file
  cat("Processing raw data...", file=log_file_name, append=TRUE, sep="\n")
  cat("Processing raw data...", sep="\n")
  data=Clean_Up_Input(rawData, groups, log_file_name=log_file_name)
  og_processed_data=data$processed
  groups=data$groups

  #Centroids or medoids?
  if(centroids==TRUE){
    centroid_flag=TRUE
  }else{
    centroid_flag=FALSE
  }

  #Original data heatmap
  if(heatmap==TRUE){
    cat("Creating original data heatmap...", file=log_file_name, append=TRUE, sep="\n")
    cat("Creating original data heatmap...", sep="\n")
    breaks=seq(0, #start point of color key
               as.numeric(quantile(data.matrix(data$processed[2:nrow(data$processed), 2:ncol(data$processed)]), 0.99)),  #end point of color key
               by=0.05) #length of sub-division
    mycol <- colorpanel(n=length(breaks)-1, low="black", high= "yellow") #heatmap colors
    suppressWarnings(DDheatmap(data.matrix(data$processed[2:nrow(data$processed), 2:ncol(data$processed)]), #the data matrix
                               Colv=FALSE, # No clustering of columns
                               Rowv = FALSE, #no clustering of rows
                               dendrogram = "none", #do not generate dendrogram
                               col=mycol, #colors used in heatmap
                               ColSideColors = as.color(Renumber(data$processed[1,2:ncol(data$processed)]), alpha=1, seed=4), #column color bar
                               RowSideColors = as.color(Renumber(data$processed[2:nrow(data$processed),1]), alpha=1, seed=2), # row color bar
                               breaks=breaks, #color key details
                               trace="none", #no trace on map
                               na.rm=TRUE, #ignore missing values
                               margins = c(5,5), # size and layout of heatmap window
                               labRow=NA, #turn off gene labels
                               labCol=NA, #turn off cell labels
                               xlab = "Samples", #x axis title
                               ylab =  "Genes", # y axis title
                               main = paste0("Original data: ", filename))) #main title
    }

  #Remove cell cycle gene cluster (optional)
  if(removeCC==TRUE){
    cat("Removing cell cycle clusters...", file=log_file_name, append=TRUE, sep="\n")
    cat("Removing cell cycle clusters...", sep="\n")
    data=Remove_Cell_Cycle(data$processed, species, log_file_name)
  }else{
    data=data$processed
  }
  if(write==TRUE){
    write.table(data, paste0(location, "data_processed_", filename, ".txt"), sep="\t")
    write.table(groups, paste0(location, "groups_processed_", filename, ".txt"), sep="\t")
  }

  #Calculate medoids, medoid correlations, blacklist to create new combine medoids
  cat("Combining similar clusters...", file=log_file_name, append=TRUE, sep="\n")
  cat("Combining similar clusters...", sep="\n")
  BL=Blacklist_Groups(data, groups, rhop, centroid_flag, log_file_name)
  newMedoids=BL$newMedoids
  groupsMedoids=BL$newGroups

  #Create synthetic doublets to get average synthetic profiles
  cat("Creating synthetic doublet profiles...", file=log_file_name, append=TRUE, sep="\n")
  cat("Creating synthetic doublet profiles...", sep="\n")
  if(.Platform$OS.type=="unix"){
    sink("/dev/null") #hides DeconRNASeq output
    synthProfilesx=Synthetic_Doublets(data, groups, groupsMedoids, newMedoids, num_doubs, log_file_name=log_file_name, only50=only50, location=location)
    sink()
  }else{
    synthProfilesx=Synthetic_Doublets(data, groups, groupsMedoids, newMedoids, num_doubs, log_file_name=log_file_name, only50=only50, location=location)
  }
  synthProfiles=synthProfilesx$averagesAverages
  doubletCellsInput2=synthProfilesx$doubletCellsInput2
  if(write==TRUE){
    write.table(doubletCellsInput2, paste0(location, "Synth_doublet_info_", filename, ".txt"), sep="\t")
  }

  #Calculate doublets using DeconRNASeq
  cat("Step 1: Removing possible doublets...", file=log_file_name, append=TRUE, sep="\n")
  cat("Step 1: Removing possible doublets...", sep="\n")
  print('data')
  print(data)
  lapply(names(data), function(d) { data[[d]][!is.finite(data[[d]])] <-0; data[[d]] })
  print('newMedoids')
  print(newMedoids)
  if(.Platform$OS.type=="unix"){
    sink("/dev/null") #hides DeconRNASeq output
    doubletTable=Is_A_Doublet(data, newMedoids, groups, synthProfiles, log_file_name=log_file_name)
    sink()
  }else{

    doubletTable=Is_A_Doublet(data, newMedoids, groups, synthProfiles, log_file_name=log_file_name)
  }
  if(write==TRUE){
    write.table(doubletTable$isADoublet, paste0(location, "DRS_doublet_table_", filename, ".txt"), sep="\t")
    write.table(doubletTable$resultsreadable, paste0(location, "DRS_results_", filename, ".txt"), sep="\t")
  }

  #Recluster doublets and non-doublets
  cat("Step 2: Re-clustering possible doublets...", file=log_file_name, append=TRUE, sep="\n")
  cat("Step 2: Re-clustering possible doublets...", sep="\n")
  reclusteredData=Recluster(isADoublet=doubletTable$isADoublet, data, groups, log_file_name = log_file_name)
  data=reclusteredData$newData2$processed
  groups=reclusteredData$newData2$groups
  write.table(data, paste0(location, "data_processed_reclust_", filename, ".txt"), sep="\t", col.names = NA, quote=FALSE)
  write.table(groups, paste0(location, "groups_processed_reclust_", filename, ".txt"), sep="\t")

  #Run Pseudo Marker Finder to identify clusters with no unique gene expression
  if(PMF==FALSE){
    cat("SKIPPING Step 3: Rescuing cells with unique gene expression...", file=log_file_name, append=TRUE, sep="\n")
    cat("SKIPPING Step 3: Rescuing cells with unique gene expression...", sep="\n")
    PMFresults=NULL
  }else{
    cat("Step 3: Rescuing cells with unique gene expression...", file=log_file_name, append=TRUE, sep="\n")
    cat("Step 3: Rescuing cells with unique gene expression...", sep="\n")
    if(useFull==TRUE){
      PMFresults=Pseudo_Marker_Finder(as.data.frame(groups), redu_data2=paste0(location, "data_processed_reclust_", filename, ".txt"), full_data2=fullDataFile, min_uniq=min_uniq, log_file_name=log_file_name, nCores=nCores)
    }else{
      PMFresults=Pseudo_Marker_Finder(as.data.frame(groups), redu_data2=paste0(location, "data_processed_reclust_", filename, ".txt"), full_data2=NULL, min_uniq=min_uniq, log_file_name=log_file_name, nCores=nCores)
    }
    if(write==TRUE){
      write.table(PMFresults, paste0(location, "new_PMF_results_", filename, ".txt"), sep="\t")
    }
  }


  #Doublet Detection method 2: Pseudo_Marker_Finder
  allClusters=unique(groups[,1])
  if(PMF==FALSE){
    newDoubletClusters=allClusters
  }else{
    hallmarkClusters=as.numeric(unique(PMFresults[,2]))
    newDoubletClusters=setdiff(allClusters, hallmarkClusters)
  }


  #Doublet Detection method 1: Is_A_Doublet
  uniqueClusters=as.character(unique(groups[,2]))
  DeconCalledFreq=as.data.frame(matrix(nrow=length(allClusters), ncol=1), row.names = uniqueClusters)
  for(clus in 1:length(allClusters)){ #modified this line, was originally "clus in allClusters"
    temp1=subset(doubletTable$isADoublet, Group_Cluster==uniqueClusters[clus])
    if(nrow(temp1)==0){ #not an original cluster, only a new doublet cluster
      DeconCalledFreq[clus,1]=100
    }else{
      DeconCalledFreq[clus,1]=(length(which(temp1$isADoublet==TRUE))/nrow(temp1))*100
    }
  }

  #Combine to find real doublets
  if(PMF==FALSE){
    finalDoublets=row.names(doubletTable$isADoublet)[which(doubletTable$isADoublet$isADoublet==TRUE)] #this gives you the names of cells called as doublets by deconvolution
  }else{
    finalDoublets=intersect(row.names(doubletTable$isADoublet)[which(doubletTable$isADoublet$isADoublet==TRUE)],row.names(subset(groups, groups[,1] %in% newDoubletClusters)))
  }

  #Results
  finalDoubletCellCall=groups[row.names(groups) %in% finalDoublets,]
  finalNotDoubletCellCall=groups[!(row.names(groups) %in% finalDoublets),]
  if(write==TRUE){
    write.table(finalDoubletCellCall, paste0(location, "Final_doublets_groups_", filename, ".txt"), sep="\t")
    write.table(finalNotDoubletCellCall, paste0(location, "Final_nondoublets_groups_", filename, ".txt"), sep="\t")
  }

  #Subset expression matrix for doublets and save
  doublets_matrix=cbind(og_processed_data[,1],og_processed_data[,which(colnames(og_processed_data) %in% row.names(finalDoubletCellCall))])
  if(write==TRUE){
    write.table(doublets_matrix, paste0(location, "Final_doublets_exp_", filename, ".txt"), sep="\t")
  }

  #Heatmap of cells removed as doubets
  if(heatmap==TRUE){
    cat("Creating doublets heatmap...", file=log_file_name, append=TRUE, sep="\n")
    cat("Creating doublets heatmap...", sep="\n")
    breaks=seq(0, #start point of color key
               as.numeric(quantile(data.matrix(doublets_matrix[2:nrow(doublets_matrix), 2:ncol(doublets_matrix)]), 0.99)),  #end point of color key
               by=0.05) #length of sub-division
    mycol <- colorpanel(n=length(breaks)-1, low="black", high= "yellow") #heatmap colors
    suppressWarnings(DDheatmap(data.matrix(doublets_matrix[2:nrow(doublets_matrix), 2:ncol(doublets_matrix)]), #the data matrix
                               Colv=FALSE, # No clustering of columns
                               Rowv = FALSE, #no clustering of rows
                               col=mycol, #colors used in heatmap
                               dendrogram="none", #turn of dendrogram generation
                               ColSideColors = as.color(Renumber(doublets_matrix[1,2:ncol(doublets_matrix)]), alpha=1, seed=4), #column color bar
                               RowSideColors = as.color(Renumber(doublets_matrix[2:nrow(doublets_matrix),1]), alpha=1, seed=2), # row color bar
                               breaks=breaks, #color key details
                               trace="none", #no trace on map
                               na.rm=TRUE, #ignore missing values
                               margins = c(5,5), # size and layout of heatmap window
                               labRow=NA, #turn off gene labels
                               labCol=NA, #turn off cell labels
                               xlab = "Samples", #x axis title
                               ylab =  "Genes", # y axis title
                               main = paste0("Doublets: ", filename))) #main title)
  }

  #Subset expression matrix for non-doublets and save
  nondoublets_matrix=cbind(og_processed_data[,1],og_processed_data[,which(colnames(og_processed_data) %in% row.names(finalNotDoubletCellCall))])
  if(write==TRUE){
    write.table(nondoublets_matrix, paste0(location, "Final_nondoublets_exp_", filename, ".txt"), sep="\t")
  }

  #New heatmap of non-doublet cells
  if(heatmap==TRUE){
    cat("Creating non-doublets heatmap...", file=log_file_name, append=TRUE, sep="\n")
    cat("Creating non-doublets heatmap...", sep="\n")
    breaks=seq(0, #start point of color key
               as.numeric(quantile(data.matrix(nondoublets_matrix[2:nrow(nondoublets_matrix), 2:ncol(nondoublets_matrix)]), 0.99)),  #end point of color key
               by=0.05) #length of sub-division
    mycol <- colorpanel(n=length(breaks)-1, low="black", high= "yellow") #heatmap colors
    suppressWarnings(DDheatmap(data.matrix(nondoublets_matrix[2:nrow(nondoublets_matrix), 2:ncol(nondoublets_matrix)]), #the data matrix
                               Colv=FALSE, # No clustering of columns
                               Rowv = FALSE, #no clustering of rows
                               col=mycol, #colors used in heatmap
                               dendrogram="none", #turn of dendrogram generation
                               ColSideColors = as.color(Renumber(nondoublets_matrix[1,2:ncol(nondoublets_matrix)]), alpha=1, seed=4), #column color bar
                               RowSideColors = as.color(Renumber(nondoublets_matrix[2:nrow(nondoublets_matrix),1]), alpha=1, seed=2), # row color bar
                               breaks=breaks, #color key details
                               trace="none", #no trace on map
                               na.rm=TRUE, #ignore missing values
                               margins = c(5,5), # size and layout of heatmap window
                               labRow=NA, #turn off gene labels
                               labCol=NA, #turn off cell labels
                               xlab = "Samples", #x axis title
                               ylab =  "Genes", # y axis title
                               main = paste0("Non-Doublets: ", filename))) #main title
  }

  #last message
  cat("Finished!", file=log_file_name, append=TRUE, sep="\n")
  cat("Finished!", sep="\n")
  
  #close the log file connection
  close(log_con)

  return(list(data_processed=data,
              groups_processed=groups,
              DRS_doublet_table=doubletTable$isADoublet,
              DRS_results=doubletTable$resultsreadable,
              PMF_results=PMFresults,
              Final_doublets_groups=finalDoubletCellCall,
              Final_nondoublets_groups=finalNotDoubletCellCall,
              Final_doublets_exp=doublets_matrix,
              Final_nondoublets_exp=nondoublets_matrix,
              Synth_doublet_info=doubletCellsInput2))

}

## Run Doublet Decon ##
results <- Main_Doublet_Decon(rawDataFile = processed$newExpressionFile, 
  groupsFile = processed$newGroupsFile, 
  filename = "DoubletDecon_results",
  location = paste0(args$out, "/"),
  fullDataFile = NULL, 
  removeCC = args$removeCC, 
  species = args$species, 
  rhop = args$rhop,
  write = TRUE, 
  PMF = args$pmf, 
  useFull = FALSE, 
  heatmap = args$heatmap, 
  centroids=args$centroids, 
  num_doubs=args$num_doubs, 
  only50=args$only50, 
  min_uniq=args$min_uniq, 
  nCores = args$nCores)




doublets <- read.table(paste0(args$out, "/Final_doublets_groups_DoubletDecon_results.txt"))
doublets$Barcode <- gsub("\\.", "-",rownames(doublets))
doublets$DoubletDecon_DropletType <- "doublet"
doublets$V1 <- NULL
doublets$V2 <- NULL


singlets <- read.table(paste0(args$out, "/Final_nondoublets_groups_DoubletDecon_results.txt"))
singlets$Barcode <- gsub("\\.", "-",rownames(singlets))
singlets$DoubletDecon_DropletType <- "singlet"
singlets$V1 <- NULL
singlets$V2 <- NULL

doublets_singlets <- rbind(singlets,doublets)

fwrite(doublets_singlets, paste0(args$out, "/DoubletDecon_doublets_singlets.tsv"), sep = "\t", append = FALSE)


### Make a summaruy of the number of singlets and doublets
summary <- as.data.frame(table(doublets_singlets$DoubletDecon_DropletType))
colnames(summary) <- c("Classification", "Droplet N")
fwrite(summary, paste0(args$out,"/DoubletDecon_doublet_summary.tsv"), sep = "\t", append = FALSE)