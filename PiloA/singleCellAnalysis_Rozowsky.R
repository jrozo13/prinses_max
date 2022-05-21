# Single-cell Pilocytic Astrocytoma analysis (Alex's Lemonade Stand)
# Last updated: 26/04/2022

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
load(paste0(wd, "PA/PA_Data/Rozowsky.SingleCellExpObject_19.05.2022.RData"))
library(Seurat)
library(SingleCellExperiment)
Rozowsky_seurat <- as.Seurat(Rozowsky_sce)

object.list <- SplitObject(Rozowsky_seurat, split.by = "orig.ident")
for (i in 1:length(object.list)) {
  print(paste0("index: ", i, "; nCells: ", length(Cells(object.list[[i]]))))
}
# remove samples with < 30 cells after filtering

for (i in 1:length(object.list)) {
  object.list[[i]] <- NormalizeData(object.list[[i]], verbose = F)
  object.list[[i]] <- FindVariableFeatures(object.list[[i]], nfeatures = 4000, verbose = F)
}

# Find anchors
anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:30, verbose = F)
# Integrate data
integrated <- IntegrateData(anchorset = anchors, verbose = F)

integrated <- ScaleData(integrated, verbose = F)
integrated <- RunPCA(integrated, features = VariableFeatures(integrated))
ElbowPlot(object = integrated, ndims = 50)
# use 30 PCs

dims.use <- 30
integrated <- RunTSNE(integrated, dims = 1:dims.use, verbose = F)
integrated <- RunUMAP(integrated, dims = 1:dims.use, verbose = F)
# assign(paste0("integrated_mito", mito_fraction.check), integrated)
DimPlot(integrated, reduction = "umap", group.by = "orig.ident")
DimPlot(integrated, reduction = "umap", group.by = "updated_location")

ST_genes <- c("FOXG1", "SIX6", "VAX1", "SIX3", "C14orf39")
PF_genes <-
FeaturePlot(integrated,
            reduction = "umap",
            split.by = "updated_location",
            features = ST_genes)

# We do not see separation by patient or batches
FeaturePlot(integrated,
            reduction = "tsne", 
            features = c("nUMI", "nGene", "percent.mito"))

########## Clustering ##########
integrated <- FindNeighbors(integrated,
                           dims = 1:dims.use,
                           annoy.metric = "euclidean")
integrated <- FindClusters(integrated, 
                           resolution = 0.2, 
                           group.singletons = FALSE,
                           algorithm = 1)
DimPlot(integrated,
        pt.size = 1,
        reduction = "umap")

library(SingleR)
library(celldex)
ref.data <- celldex::HumanPrimaryCellAtlasData(ensembl = FALSE)
# ensembl == TRUE: will use the ENSEMBL IDs
cluster_predictions <- SingleR(as.SingleCellExperiment(integrated), 
                               ref = ref.data,
                               labels = ref.data$label.main,
                               clusters = integrated$integrated_snn_res.0.2)
cluster_idents <- data.frame(cluster_predictions@rownames, cluster_predictions@listData$labels)
cell_idents <- merge(x = integrated@meta.data %>% rownames_to_column() %>%  select(rowname, integrated_snn_res.0.2), 
                     y = cluster_idents,
                     by.x = "integrated_snn_res.0.2",
                     by.y = "cluster_predictions.rownames")%>%
  column_to_rownames(var = "rowname")
cell_idents <- cell_idents[rownames(integrated@meta.data),]
integrated <- AddMetaData(integrated, metadata = cell_idents$cluster_predictions.listData.labels, col.name = "SingleR.Annotation")
DimPlot(integrated, reduction = "umap", group.by = "SingleR.Annotation", pt.size = 1.5)
DimPlot(integrated, reduction = "umap", group.by = "integrated_snn_res.0.2")

markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

DoHeatmap(integrated, features = top5$gene) + NoLegend()

# Immune cells vs glioma vs stromal cells
FeaturePlot(integrated, features = "PTPRC")

assign("Rozowsky_seurat", integrated)
save(Rozowsky_seurat, file = paste0(wd, "PA/PA_Data/Rozowsky.SCEAnalaysis_19.05.2022.RData"))


