# Single-cell glioma microenvironment analysis (Pombo-Antunes et al, 2021)
# Last updated: 17/03/2022

########## Initialize ##########
## Install global packages ##
library(dplyr)
library(tidyr)
library(tidyverse)
library(readxl)
library(RColorBrewer)
library(circlize)
col_annotations <- list(
  Sex = c("M" = "skyblue", "F" = "pink"),
  Location = c("Spinal" = "#EF553B", 
               "Posterior fossa" = "#636EFA", 
               "Supratentorial" = "#00CC96"),
  Age = colorRamp2(c(0,20), c("#edf8e9", "#74c476"))
)

wd <- "/Users/jrozowsky/Library/Mobile Documents/com~apple~CloudDocs/Documents/PMC/"
setwd(wd)
fwd <- paste0(wd, "PA/Analysis/Figures/") # figure working directory

########## Load Pombo dataset ##########
library(Seurat)
library(SingleCellExperiment)

# Compile all patient IDs. Make Seurat objects with for-loop
cell_metaData <- read.csv(file = "Data/Pombo_2021/annot_Human_ND_GBM_Full.csv")
file_dir <- paste0(wd, "Data/Pombo_2021/filtered_feature_bc_matrix_HumanNewlyDiagnGBM/filtered_feature_bc_matrix/")
counts <- Read10X(data.dir = file_dir)
pombo_all <- CreateSeuratObject(counts = counts, assay = "RNA")

cell_metaData$cell %in% Cells(pombo_all) %>% table()
# all newly diagnosed cells are in seurat object -> subset for those cells
pombo_all$cellID <- Cells(pombo_all)
pombo <- subset(x = pombo_all, subset = cellID %in% cell_metaData$cell)
rm(pombo_all)

pombo@meta.data <- cbind(pombo@meta.data, cell_metaData %>% column_to_rownames(var = "cell"))
identical(rownames(pombo@meta.data), pombo@meta.data$cellID)
# check that the merged data frames were in order

########## Quality Control ########## 
library(scater)
save(pombo, file = paste0(wd, "PA/PA_Data/Pombo.SeuratObject_17.03.2022.RData"))
mito.genes <- grep("^MT-", rownames(pombo@assays$RNA@counts))
pombo$percent.mito <- Matrix::colSums(pombo@assays$RNA@counts[mito.genes, ])/Matrix::colSums(pombo@assays$RNA@counts)

pombo$nCount.outlier.low <- isOutlier(pombo$nCount_RNA, nmads=3, type="lower", log=TRUE)
pombo$nFeature.outlier.low <- isOutlier(pombo$nFeature_RNA, nmads=3, type="lower", log=TRUE)
pombo$mito.outlier.high <- isOutlier(pombo$percent.mito, nmads=3, type="higher", log=TRUE)
# dont remove any cells

is.mito <- grepl("^MT-", rownames(pombo))
sce_clean <- addPerCellQC(pombo@assays$RNA@counts, subsets=list(Mt=is.mito), flatten=T)

pombo <- NormalizeData(pombo,verbose = F)
pombo <- FindVariableFeatures(pombo,verbose=F)
pombo <- ScaleData(pombo,verbose=F)
pombo <- RunPCA(pombo, features =VariableFeatures(pombo))
ElbowPlot(object = pombo, ndims = 50)
pombo <- FindNeighbors(pombo, dims = 1:30, annoy.metric = "euclidean")
pombo <- FindClusters(pombo, resolution = 0.1, group.singletons = TRUE)
pombo <- RunUMAP(pombo, dims = 1:30, verbose=F)
DimPlot(object = pombo, group.by = "seurat_clusters")
DimPlot(object = pombo, group.by = "seurat_clusters")

########## Load Pombo dataset ##########
library(SingleR)
ref.pombo <- as.SingleCellExperiment(pombo, assay = "RNA")

