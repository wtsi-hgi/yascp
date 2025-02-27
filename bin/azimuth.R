#!/usr/bin/env Rscript

## adapted from https://azimuth.hubmapconsortium.org

## expects input *.h5 file as argument
## writes plots to a single PDF file Rplots.pdf


#!/usr/bin/env Rscript

## adapted from https://azimuth.hubmapconsortium.org

## expects input *.h5 file as argument
## writes plots to a single PDF file Rplots.pdf
library(Azimuth)
library(Seurat)
library(SeuratDisk)
library(Matrix)
library(hdf5r)
library(ggplot2)
library(tools)
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
# reference files are expected in the following directory

NNTransform <- function(
  object,
  meta.data,
  neighbor.slot = "query_ref.nn",
  key = 'ori.index'
) {
  on.exit(expr = gc(verbose = FALSE))
  ind <- Indices(object[[neighbor.slot]])
  ori.index <- t(x = sapply(
    X = 1:nrow(x = ind),
    FUN = function(i) {
      return(meta.data[ind[i, ], key])
    }
  ))
  rownames(x = ori.index) <- rownames(x = ind)
  slot(object = object[[neighbor.slot]], name = "nn.idx") <- ori.index
  return(object)
}

#######################
######## main #########
#######################

# inputfile.h5ad = './AZ_pre_QC_adata_full.h5ad'

cat("inputfile.h5ad = ", inputfile.h5ad, "\n")
# input_dir = args[1]

options("Azimuth.map.ndims" = 50) # used in helpers.R::LoadReference()

# Ensure Seurat v4.0 or higher is installed
if (packageVersion(pkg = "Seurat") < package_version(x = "4.0.0")) {
  stop("Mapping datasets requires Seurat v4 or higher.", call. = FALSE)
}

# Ensure glmGamPoi is installed
if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    BiocManager::install("glmGamPoi")
  }
}

# copy input file and convert to h5seurat
Convert(inputfile.h5ad, dest="h5seurat", overwrite = TRUE)
inputfile.h5seurat <- paste0(file_path_sans_ext(inputfile.h5ad), ".h5seurat")
cat("inputfile.h5seurat = ", inputfile.h5seurat, "\n")
cat("Loading file", inputfile.h5seurat, "\n")

query <- LoadH5Seurat(inputfile.h5seurat, meta.data = FALSE, misc = FALSE)
# Download the Azimuth reference and extract the archive
saveRDS(query, file = "query.rds")
# Load the reference
# Change the file path based on where the reference is located on your system.
## reference <- LoadReference(path = "https://seurat.nygenome.org/azimuth/references/v1.0.0/human_pbmc")
reference <- Azimuth:::LoadReference(REFERENCE_DIR)

cat("query file loaded.\n")


# Calculate percent mitochondrial genes if the query contains genes
# matching the regular expression "^MT-"
if (any(grepl(pattern = '^MT-', x = rownames(x = query)))) {
  query <- PercentageFeatureSet(
    object = query,
    pattern = '^MT-',
    col.name = 'percent.mt',
    assay = "RNA"
  )
}

# Preprocess with SCTransform
query <- SCTransform(
  object = query,
  assay = "RNA",
  new.assay.name = "refAssay",
  residual.features = rownames(x = reference$map),
  reference.SCT.model = reference$map[["refAssay"]]@SCTModel.list$refmodel,
  method = 'glmGamPoi',
  ncells = 2000,
  n_genes = 2000,
  do.correct.umi = FALSE,
  do.scale = FALSE,
  do.center = TRUE
)

# Find anchors between query and reference
anchors <- FindTransferAnchors(
  reference = reference$map,
  query = query,
  k.filter = NA,
  reference.neighbors = "refdr.annoy.neighbors",
  reference.assay = "refAssay",
  query.assay = "refAssay",
  reference.reduction = "refDR",
  normalization.method = "SCT",
  features = intersect(rownames(x = reference$map), VariableFeatures(object = query)),
  dims = 1:50,
  n.trees = 20,
  mapping.score.k = 100
)

# Transfer cell type labels and impute protein expression
#
# Transferred labels are in metadata columns named "predicted.*"
# The maximum prediction score is in a metadata column named "predicted.*.score"
# The prediction scores for each class are in an assay named "prediction.score.*"
# The imputed assay is named "impADT" if computed

# levels<- list("celltype.l2", "celltype.l1", "celltype.l3")
for (celltype_level in levels) {
    id <- celltype_level  # Choose a correct metadata column (class, subclass, etc.)
    print(paste("Processing:", id))

    # Ensure the selected column exists
    if (!id %in% colnames(reference$map[[]])) {
        stop(paste("Column", id, "not found in reference metadata. Available columns:", 
                   paste(colnames(reference$map[[]]), collapse = ", ")))
    }

    # Extract correct cell type labels
      print(id)
      refdata <- lapply(X = celltype_level, function(x) {
        reference$map[[x, drop = TRUE]]
      })

      names(x = refdata) <- celltype_level
      if (FALSE) {
        refdata[["impADT"]] <- GetAssayData(
          object = reference$map[['ADT']],
          slot = 'data'
        )
      }
      query <- TransferData(
        reference = reference$map,
        query = query,
        dims = 1:50,
        anchorset = anchors,
        refdata = refdata,
        n.trees = 20,
        store.weights = TRUE
      )

      # Calculate the embeddings of the query data on the reference SPCA
      query <- IntegrateEmbeddings(
        anchorset = anchors,
        reference = reference$map,
        query = query,
        reductions = "pcaproject",
        reuse.weights.matrix = TRUE
      )

      # Calculate the query neighbors in the reference
      # with respect to the integrated embeddings
      query[["query_ref.nn"]] <- FindNeighbors(
        object = Embeddings(reference$map[["refDR"]]),
        query = Embeddings(query[["integrated_dr"]]),
        return.neighbor = TRUE,
        l2.norm = TRUE
      )

      # The reference used in the app is downsampled compared to the reference on which
      # the UMAP model was computed. This step, using the helper function NNTransform,
      # corrects the Neighbors to account for the downsampling.
      query <- NNTransform(
        object = query,
        meta.data = reference$map[[]]
      )

      # Project the query to the reference UMAP.
      query[["proj.umap"]] <- RunUMAP(
        object = query[["query_ref.nn"]],
        reduction.model = reference$map[["refUMAP"]],
        reduction.key = 'UMAP_'
      )


      # Calculate mapping score and add to metadata
      query <- AddMetaData(
        object = query,
        metadata = MappingScore(anchors = anchors),
        col.name = paste0("mapping.score.", id)
      )


      # VISUALIZATIONS
      cat("Generating Visualizations ...\n")
      pdf(NULL) # switch off automatic generation of Rplots.pdf
      predicted.id <- paste0("predicted.", id)
      predicted.id.score <- paste0(predicted.id, ".score")

      # write a table of cell-type assignments, prediction and mapping scores:
      fnam.table <- paste0(prefix,"___",gsub(".", "_", predicted.id, fixed = TRUE),".tsv")
      data <- FetchData(object = query, vars = c(predicted.id, predicted.id.score, paste0("mapping.score.", id)), slot = "data")
      write.table(data, fnam.table, quote = FALSE, sep="\t")
      #gzip(fnam.table, overwrite = TRUE)

      # make a barplot for the number of cells in each celltype category
      ctyp.counts <- function(vec, label) {
        tb <- table(as.factor(vec))
        v <- vector(mode = "character", length = length(tb))
        v[] <- label
        tabf <- cbind(v, as.data.frame(tb))
        names(tabf) <- c("threshold", "cell_type", "count")
        tabf
      }
      tdf <- ctyp.counts(data[,1], label = "all")
      for (thresh in c(.5, .8)) {
        tag <- paste0("score>",format(thresh, digits=1))
        tdf <- rbind(
          tdf,
          ctyp.counts(data[data[,2] > thresh, 1],
            label = paste0("score>",format(thresh, digits=1))
          )
        )
      }
      p <- ggplot(tdf, aes(x=cell_type, y=count, fill = threshold))
      p <- p + geom_col(position = "dodge")
      p <- p + theme(axis.text.x = element_text(angle = 90))
      ggsave(paste0(prefix,"___",id,".ncells_by_type_barplot.pdf"))

      # DimPlot of the reference
      #ref.plt <- DimPlot(object = reference$plot, reduction = "refUMAP", group.by = id, label = TRUE) + NoLegend()
      #ggsave(file = "ref_umap.pdf", plot = ref.plt)
      #ggsave("ref_umap.pdf")

      # DimPlot of the query, colored by predicted cell type
      DimPlot(object = query, reduction = "proj.umap", group.by = predicted.id, label = TRUE) + NoLegend()
      ggsave(paste0(prefix,"___",id,".query_umap.pdf"))

      # Plot the score for the predicted cell type of the query
      FeaturePlot(object = query, features = paste0(predicted.id, ".score"), reduction = "proj.umap")
      ggsave(paste0(prefix,"___",id,".prediction_score_umap.pdf"))
      VlnPlot(object = query, features = paste0(predicted.id, ".score"), group.by = predicted.id) + NoLegend()
      ggsave(paste0(prefix,"___",id,".prediction_score_vln.pdf"))

      # Plot the mapping score
      FeaturePlot(object = query, features = paste0("mapping.score.", id), reduction = "proj.umap")
      ggsave(paste0(prefix,"___",id,".mapping_score_umap.pdf"))
      VlnPlot(object = query, features = paste0("mapping.score.", id), group.by = predicted.id) + NoLegend()
      ggsave(paste0(prefix,"___",id,".mapping_score_vln.pdf"))
    # Continue processing...
}


  # save mapped data set
  #save(query, file = "azimuth.bin")
  # saveRDS(query, file = paste0(prefix,"azimuth.rds"))
  # load("azimuth.bin" )