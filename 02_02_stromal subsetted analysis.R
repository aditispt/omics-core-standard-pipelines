###############################################################
# Project: Sex-specific decidual stromal cell analysis
#
# Script: 05_Recluster_Stromal.R
#
# Purpose:
# Independently recluster stromal cells from each group to
# identify stromal subpopulations.
###############################################################

##############################
## Load libraries
##############################

library(Seurat)
library(ggplot2)
library(patchwork)

##############################
## Load stromal objects
##############################

HF_stroma <- readRDS("Seurat_Objects/HF_stroma.rds")
HM_stroma <- readRDS("Seurat_Objects/HM_stroma.rds")
MF_stroma <- readRDS("Seurat_Objects/MF_stroma.rds")
MM_stroma <- readRDS("Seurat_Objects/MM_stroma.rds")

##############################
## Create folders
##############################

dir.create("Figures/Stroma_Recluster", recursive = TRUE, showWarnings = FALSE)
dir.create("Results/Stroma_Recluster", recursive = TRUE, showWarnings = FALSE)

###############################################################
## Function to recluster stromal cells
###############################################################

process_stroma <- function(seurat_object, sample_name){
  
  cat("\n=====================================\n")
  cat("Processing:", sample_name, "\n")
  cat("Cells:", ncol(seurat_object), "\n")
  cat("=====================================\n")
  
  #############################################################
  ## RNA preprocessing
  #############################################################
  
  DefaultAssay(seurat_object) <- "RNA"
  
  seurat_object <- NormalizeData(
    seurat_object,
    normalization.method = "LogNormalize",
    scale.factor = 10000
  )
  
  seurat_object <- FindVariableFeatures(
    seurat_object,
    selection.method = "vst",
    nfeatures = 2000
  )
  
  #############################################################
  ## Scaling
  #############################################################
  
  seurat_object <- ScaleData(
    seurat_object,
    features = rownames(seurat_object)
  )
  
  #############################################################
  ## PCA
  #############################################################
  
  seurat_object <- RunPCA(
    seurat_object,
    features = VariableFeatures(seurat_object)
  )
  
  #############################################################
  ## Save PCA plots
  #############################################################
  
  pdf(
    paste0("Figures/Stroma_Recluster/",
           sample_name,
           "_PCA.pdf"),
    width = 8,
    height = 6
  )
  
  print(
    DimPlot(
      seurat_object,
      reduction = "pca"
    ) +
      ggtitle(sample_name)
  )
  
  dev.off()
  
  #############################################################
  ## Elbow plot
  #############################################################
  
  pdf(
    paste0("Figures/Stroma_Recluster/",
           sample_name,
           "_Elbow.pdf"),
    width = 6,
    height = 5
  )
  
  print(
    ElbowPlot(
      seurat_object,
      ndims = 30
    )
  )
  
  dev.off()
  
  #############################################################
  ## Clustering
  #############################################################
  
  seurat_object <- FindNeighbors(
    seurat_object,
    dims = 1:15
  )
  
  seurat_object <- FindClusters(
    seurat_object,
    resolution = 0.3
  )
  
  #############################################################
  ## UMAP
  #############################################################
  
  seurat_object <- RunUMAP(
    seurat_object,
    dims = 1:15
  )
  
  #############################################################
  ## Save UMAP
  #############################################################
  
  pdf(
    paste0("Figures/Stroma_Recluster/",
           sample_name,
           "_UMAP.pdf"),
    width = 8,
    height = 6
  )
  
  print(
    DimPlot(
      seurat_object,
      reduction = "umap",
      label = TRUE,
      repel = TRUE
    ) +
      ggtitle(sample_name)
  )
  
  dev.off()
  
  #############################################################
  ## Find markers
  #############################################################
  
  markers <- FindAllMarkers(
    seurat_object,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25
  )
  
  write.csv(
    markers,
    paste0(
      "Results/Stroma_Recluster/",
      sample_name,
      "_ClusterMarkers.csv"
    ),
    row.names = FALSE
  )
  
  #############################################################
  ## Save object
  #############################################################
  
  saveRDS(
    seurat_object,
    paste0(
      "Results/Stroma_Recluster/",
      sample_name,
      "_Reclustered.rds"
    )
  )
  
  return(seurat_object)
  
}

###############################################################
## Run analysis
###############################################################

HF_stroma <- process_stroma(
  HF_stroma,
  "Healthy_Female"
)

HM_stroma <- process_stroma(
  HM_stroma,
  "Healthy_Male"
)

MF_stroma <- process_stroma(
  MF_stroma,
  "Miscarriage_Female"
)

MM_stroma <- process_stroma(
  MM_stroma,
  "Miscarriage_Male"
)

###############################################################
## Finished
###############################################################

cat("\n=====================================\n")
cat("All stromal datasets processed.\n")
cat("=====================================\n")