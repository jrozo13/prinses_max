# Single-cell Pilocytic Astrocytoma analysis: Integrated data-sets
# Last updated: 19/05/2022

########## Initialize ##########
## Install global packages ##
library(dplyr)
library(tidyr)
library(tidyverse)
library(readxl)
library(RColorBrewer)
library(circlize)
library(ggplot2)
library(ggrepel)
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

########## Make Seurat Object ##########
load("PA/PA_Data/Taylor.SingleCellExpObject_19.05.2022.RData")
load("PA/PA_Data/ALS.SingleCellExpObject_19.05.2022.RData")
load("PA/PA_Data/Rozowsky.SingleCellExpObject_19.05.2022.RData")
library(Seurat)
library(SingleCellExperiment)

ALS_seurat <- as.Seurat(ALS_sce)
Taylor_seurat <- as.Seurat(Taylor_sce)
Rozowsky_seurat <- as.Seurat(Rozowsky_sce)

object.list <- c(SplitObject(ALS_seurat, split.by = "orig.ident"), 
                 SplitObject(Taylor_seurat, split.by = "orig.ident"))
for (i in 1:length(object.list)) {
  print(paste0("index: ", i, "; nCells: ", length(Cells(object.list[[i]]))))
}
# remove samples with < 30 cells after filtering
object.list[4] <- NULL

for (i in 1:length(object.list)) {
  object.list[[i]] <- NormalizeData(object.list[[i]], verbose = F)
  object.list[[i]] <- FindVariableFeatures(object.list[[i]], nfeatures = 4000, verbose = F)
}

########## Seurat alignment ##########
# Find anchors
anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:30, verbose = F)
save(anchors, file = "/Users/jrozowsky/Desktop/anchors.RData")
# Integrate data
integrated <- IntegrateData(anchorset = anchors, verbose = F)
save(integrated, file = "/Users/jrozowsky/Desktop/integrated.RData")

########## Harmony alignment ##########
pa.combined <- merge(ALS_seurat, y = c(Taylor_seurat), project = "Pilocytic Astrocytoma")
pa.combined <- pa.combined %>% 
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData(verbose = FALSE) %>% 
  RunPCA(pc.genes = pa.combined@var.genes, npcs = 20, verbose = FALSE)
pa.combined$dataset <- ifelse(grepl("SCPCS", pa.combined$ident), "ALS",
                             ifelse(grepl("PA", pa.combined$ident), "Taylor", "Rozowsky"))
DimPlot(object = pa.combined, reduction = "pca", pt.size = .1, group.by = "dataset")

pa.combined <- pa.combined %>% RunHarmony("ident", plot_convergence = TRUE)
DimPlot(object = pa.combined, reduction = "harmony", pt.size = .1, group.by = "dataset")
VlnPlot(object = pa.combined, features = "harmony_1", group.by = "dataset", pt.size = .1)

pa.combined <- pa.combined %>% 
  RunUMAP(reduction = "harmony", dims = 1:20) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:20) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

DimPlot(pa.combined, reduction = "umap", pt.size = .1)
FeaturePlot(pa.combined, features = "TMEM119", split.by = "dataset")
DimPlot(pa.combined, reduction = "umap", pt.size = .1)

library(SingleR)
library(celldex)
ref.data <- celldex::HumanPrimaryCellAtlasData(ensembl = FALSE)
# ensembl == TRUE: will use the ENSEMBL IDs
cluster_predictions <- SingleR(as.SingleCellExperiment(pa.combined), 
                               ref = ref.data,
                               labels = ref.data$label.main,
                               clusters = pa.combined$seurat_clusters)
cluster_idents <- data.frame(cluster_predictions@rownames, cluster_predictions@listData$labels)
cell_idents <- merge(x = pa.combined@meta.data %>% rownames_to_column() %>% select(rowname, seurat_clusters), 
                     y = cluster_idents,
                     by.x = "seurat_clusters",
                     by.y = "cluster_predictions.rownames") %>%
  column_to_rownames(var = "rowname")
cell_idents <- cell_idents[rownames(pa.combined@meta.data),]
pa.combined <- AddMetaData(pa.combined, metadata = cell_idents$cluster_predictions.listData.labels, col.name = "SingleR.Annotation")
DimPlot(pa.combined, reduction = "umap", group.by = "SingleR.Annotation", split.by = "updated_location")
DimPlot(pa.combined, reduction = "umap", group.by = "SingleR.Annotation", split.by = "dataset")

########## Seurat Alignment ##########
integrated$dataset <- ifelse(grepl("SCPCS", integrated$ident), "ALS",
                             ifelse(grepl("PA", integrated$ident), "Taylor", "Rozowsky"))
integrated <- ScaleData(integrated, verbose = F)
integrated <- RunPCA(integrated, features = VariableFeatures(integrated))
ElbowPlot(object = integrated, ndims = 50)
# use 20 PCs

dims.use <- 30
integrated <- RunUMAP(integrated, dims = 1:dims.use,umap.method = "uwot", verbose = F)

DimPlot(integrated, group.by = "orig.ident", pt.size = 0.05)
DimPlot(integrated, group.by = "updated_location", pt.size = 0.05)

# We do not see separation by patient or batches
FeaturePlot(integrated,
            features = c("nUMI", "nGene", "percent.mito"))
FeaturePlot(integrated,
            features = c("CD2"))

integrated <- FindNeighbors(integrated, 
                            dims = 1:dims.use,
                            annoy.metric = "euclidean")
integrated <- FindClusters(integrated, 
                           resolution = 0.5, group.singletons = TRUE,
                           algorithm = 1)
DimPlot(integrated)

library(SingleR)
library(celldex)
ref.data <- celldex::HumanPrimaryCellAtlasData(ensembl = FALSE)
# ensembl == TRUE: will use the ENSEMBL IDs
cluster_predictions <- SingleR(as.SingleCellExperiment(integrated), 
                               ref = ref.data,
                               labels = ref.data$label.main,
                               clusters = integrated$integrated_snn_res.0.5)
cluster_idents <- data.frame(cluster_predictions@rownames, cluster_predictions@listData$labels)
cell_idents <- merge(x = integrated@meta.data %>% rownames_to_column() %>%  select(rowname, integrated_snn_res.0.5), 
                     y = cluster_idents,
                     by.x = "integrated_snn_res.0.5",
                     by.y = "cluster_predictions.rownames")%>%
  column_to_rownames(var = "rowname")
cell_idents <- cell_idents[rownames(integrated@meta.data),]
integrated <- AddMetaData(integrated, metadata = cell_idents$cluster_predictions.listData.labels, col.name = "SingleR.Annotation")
DimPlot(integrated, reduction = "umap", group.by = "SingleR.Annotation")
DimPlot(integrated, reduction = "umap", group.by = "integrated_snn_res.0.5")
DimPlot(integrated, reduction = "umap", split.by = "orig.ident")
DimPlot(integrated, reduction = "umap", group.by = "dataset")


save(integrated, file = paste0(wd, "PA/PA_Data/ALS.IntegratedObject_29.04.2022.RData"))

myeloid_genes <- c("C1QA", "CD163", "C1QB", "C1QC")
t_genes <- c("CD3D", "CD3E", "CD3G", "CD247")
nk_genes <- c("KIR2DL1", "KIR2DL3", "KIR2DS4")
nksig_genes <- c("KLRC1", "GNLY", "TRDC", "FGFBP2", "KLRB1", "KLRC2")
cytotoxic_genes <- c("PRF1", "GZMB", "GZMA", "GZMH", "NKG7", "GNLY")
inhib_genes <- c("PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT")
trm_genes <- c("CA10", "ITGA1", "ITGAE", "IL2", "IL10", "CXCR6", "CXCL13", "KCNK5", "RGS1", "CRTAM", "DUSP6", "PDCD1")
cd8_genes <- c("CD8A", "KLRK1", "NKG7", "CD8B", "KLRC4", "KLRD1", "ZNF683")
cd4_genes <- c("CD4", "CD40LG", "KLRB1", "RUNX2", "IL7R")
integrated <- AddModuleScore(integrated,
                          features = list(myeloid_genes),
                          name = "ModuleScore",
                          assay = "RNA")
FeaturePlot(integrated,
            features = c("ModuleScore1"), split.by = "updated_location")

  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))
########## Annotation of immune cells ##########
# Proceed with harmony alignment
load("PA/PA_Data/Pombo.SeuratObj_30.03.2022.RData")
assign("ref.pombo", seuratObj); rm(seuratObj)
library(SingleR)
library(SingleCellExperiment)
library(Seurat)



Idents(pa.combined) <- "SingleR.Annotation"
immune.obj <- subset(x = pa.combined, idents = c("Macrophage", "Monocyte", "T_cells"))
immune.obj <- ScaleData(immune.obj, verbose = F)
immune.obj <- RunPCA(immune.obj, features = VariableFeatures(immune.obj))

immune.obj <- immune.obj %>% RunHarmony("dataset", plot_convergence = TRUE)
DimPlot(object = immune.obj, reduction = "harmony", pt.size = .1, group.by = "dataset")
VlnPlot(object = immune.obj, features = "harmony_1", group.by = "dataset", pt.size = .1)

immune.obj <- immune.obj %>%
  RunUMAP(reduction = "harmony", dims = 1:20)

save(pa.combined, integrated, file = "thisisshit.RData")
