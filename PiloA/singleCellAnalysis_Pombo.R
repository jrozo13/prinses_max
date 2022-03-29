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
sce <- SingleCellExperiment(assays=list(counts=counts))
rm(counts)

#pombo_all <- CreateSeuratObject(counts = counts, assay = "RNA")

cell_metaData$cell %in% Cells(sce) %>% table()
# all newly diagnosed cells are in seurat object -> subset for those cells
sce <- sce[, cell_metaData$cell]
sce@metadata <- data.frame(cellID = Cells(sce))
sce@metadata <- cbind(sce@metadata, cell_metaData %>% column_to_rownames(var = "cell"))
identical(rownames(sce@metadata), sce@metadata$cellID)
# check that the merged data frames were in order

########## Quality Control ########## 
library(scater)
# save(sce, file = paste0(wd, "PA/PA_Data/Pombo.SingleCellExpObject_29.03.2022.RData"))
# load(paste0(wd, "PA/PA_Data/Pombo.SingleCellExpObject_29.03.2022.RData"))
mito.genes <- grep("^MT-", rownames(sce))
sce$percent.mito <- (Matrix::colSums(counts(sce)[mito.genes, ])*100)/Matrix::colSums(counts(sce)) 
sce$nGene <- apply(counts(sce),  2,  function(x) length(x[x > 0])) # number of expressed genes
sce$nUMI <- apply(counts(sce),  2,  sum) # total UMI counts (library size)
sce$staticNr <- 1
dim(colData(sce))
colnames(colData(sce))

sce$nUMI.outlier.low <- isOutlier(sce$nUMI, nmads=3, type="lower", log=TRUE)
sum(sce$nUMI.outlier.low )

sce$nGene.outlier.low <- isOutlier(sce$nGene, nmads=3, type="lower", log=TRUE)
sum(sce$nGene.outlier.low)

sce$mito.outlier.high <- isOutlier(sce$percent.mito, nmads=3, type="higher", log=TRUE)
sum(sce$mito.outlier.high)

# Filter cells and identify outliers
metaData=as.data.frame(colData(sce))
ggplot(metaData, aes(staticNr, nUMI)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=nUMI.outlier.low)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("Total UMI counts per cell")
ggplot(metaData, aes(staticNr, nGene)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=nGene.outlier.low)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("Number of genes per cell")
ggplot(metaData, aes(staticNr, percent.mito)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=mito.outlier.high)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("% mito genes per cell")

metaData.filtered<-metaData[! (metaData$nUMI.outlier.low | metaData$nGene.outlier.low | metaData$mito.outlier.high) ,]
ggplot(metaData.filtered, aes(staticNr, nUMI)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=nUMI.outlier.low)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("Total UMI counts per cell after filtering")
ggplot(metaData.filtered, aes(staticNr, nGene)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=nGene.outlier.low)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("Number of genes per cell after filtering")
ggplot(metaData.filtered, aes(staticNr, percent.mito)) + 
  geom_violin(fill="gray80") + 
  geom_jitter(height = 0, width = 0.3, aes(col=mito.outlier.high)) +
  scale_color_manual(values=c("#00BFC4", "#F8766D"))+ggtitle("% mito genes per cell after filtering")

# Check which metaData was used to identify outliers
library(VennDiagram)
v <-venn.diagram(
  list (UMI=rownames(colData(sce)[sce$nUMI.outlier.low,]),
        gene=rownames(colData(sce)[sce$nGene.outlier.low,]),
        mito=rownames(colData(sce)[sce$mito.outlier.high,])),
  filename=NULL,
  alpha = c( 0.5,0.5,0.5),
  fill = c("green","orange","blue")
)
grid.newpage()
grid.draw(v)
rm(v)

# Remove outliers
sce_clean <- sce[,!(sce$nUMI.outlier.low | sce$nGene.outlier.low | sce$mito.outlier.high)]
dim(sce_clean)
ncol(sce_clean)-ncol(sce)

# Check how many cells were removed
sce_clean <- addPerCellQC(pombo@assays$RNA@counts, subsets=list(Mt=is.mito), flatten=T)
sce_clean <- addPerCellQC(sce_clean, flatten=T)

# Filter out by genes
ave.counts <- rowMeans(as.matrix(counts(sce_clean)))
thresh <- 0.001  #cutoff can be modified here
hist(log10(ave.counts), breaks=100, col="grey80",
     xlab = expression(Log[10]~"mean count"))
abline(v=log10(thresh), col="blue", lwd=2, lty=2)

rowData(sce_clean)$usegenes <- ave.counts>thresh
table(rowData(sce_clean)$usegenes)

# Filter out the lowly-abundant genes
sce_clean <- sce_clean[rowData(sce_clean)$usegenes, ]

########## Make Seurat Object ##########
counts <- counts(sce_clean)
rownames(counts) <- rownames(sce_clean)
colnames(counts) <- colnames(sce_clean)
seuratObj <- CreateSeuratObject(counts = counts)
meta.data.keep <- colData(sce_clean)
meta.data.keep <- meta.data.keep[,c("percent.mito")]

seuratObj <- AddMetaData(seuratObj, as.data.frame(meta.data.keep))
save(seuratObj, file = paste0(wd, "PA/PA_Data/Pombo.SeuratObj_29.03.2022.RData"))
save(sce, file = paste0(wd, "PA/PA_Data/Pombo.SingleCellExpObject_29.03.2022.RData"))

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


is.mito <-grepl("^MT-",rownames(sce_clean))
sum(is.mito)
sce_clean <- addPerCellQC(sce_clean, subsets=list(Mt=is.mito), flatten=T)

########## Load Pombo dataset ##########
library(SingleR)
ref.pombo <- as.SingleCellExperiment(pombo, assay = "RNA")

