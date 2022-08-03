# Single-cell Pilocytic Astrocytoma analysis (Rozowsky)
# Last updated: 03/08/2022

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
load(paste0(wd, "PA/PA_Data/Rozowsky.SingleCellExpObject_03.08.2022.RData"))
library(Seurat)
library(SingleCellExperiment)
sce <- as.Seurat(Rozowsky_sce); rm(Rozowsky_sce)

sce <- NormalizeData(sce)
sce <- FindVariableFeatures(sce, nfeatures = 5000)
sce <- ScaleData(sce, verbose = F)
sce <- RunPCA(sce, features = VariableFeatures(sce))
ElbowPlot(object = sce, ndims = 50)
# use 30 PCs

dims.use <- 30
sce <- RunUMAP(sce, dims = 1:dims.use, verbose = F)
DimPlot(sce, reduction = "umap", group.by = "orig.ident")
DimPlot(sce, reduction = "umap", group.by = "location")
DimPlot(sce, reduction = "umap", group.by = "date")

# see clear separation by patient --> try harmony correction

########## Harmony Alignment ##########
library(harmony)
theta.use <- 1
sce <- RunHarmony(object = sce,
                  group.by.vars = "orig.ident",
                  theta = theta.use,
                  plot_convergence = TRUE)

sce <- RunUMAP(object = sce,
               reduction = "harmony",
               dims = 1:dims.use,
               reduction.name = "umapHarmony",
               reduction.key = "umapHarmony")
DimPlot(sce, reduction = "umapHarmony", group.by = "orig.ident")
DimPlot(sce, reduction = "umapHarmony", group.by = "location")
DimPlot(sce, reduction = "umapHarmony", group.by = "date")

Rozowsky_sce.harmonyCorrection <- sce; rm(sce)
save(Rozowsky_sce.harmonyCorrection, file = paste0(wd, "PA/PA_Data/Rozowsky.SeuratHarmonyCorrection_03.08.2022.RData"))

########## Clustering ##########
load("PA/PA_Data/Rozowsky.SeuratHarmonyCorrection_03.08.2022.RData")
sce <- Rozowsky_sce.harmonyCorrection; rm(Rozowsky_sce.harmonyCorrection)

FeaturePlot(sce,
            reduction = "umapHarmony",
            split.by = "location",
            features = c("nCount_RNA", "nFeature_RNA"))

sce <- FindNeighbors(sce,
                     dims = 1:dims.use,
                     reduction = "harmony",
                     annoy.metric = "euclidean")
sce <- FindClusters(sce, 
                    resolution = 0.2, 
                    group.singletons = FALSE,
                    algorithm = 1)
DimPlot(sce,
        pt.size = 1,
        reduction = "umapHarmony",
        split.by = "location")

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
DimPlot(Rozowsky_seurat, reduction = "umap", group.by = "integrated_snn_res.0.2")

markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

DoHeatmap(integrated, features = top5$gene) + NoLegend()

# Immune cells vs glioma vs stromal cells
DimPlot(Rozowsky_seurat, group.by = "Phase")
FeaturePlot(Rozowsky_seurat, features = "H2AZ1")

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(object = Rozowsky_seurat) <- "RNA"
Rozowsky_seurat <- CellCycleScoring(Rozowsky_seurat, s.features = s.genes, g2m.features = g2m.genes,  set.ident = TRUE)
DimPlot(Rozowsky_seurat, group.by = "Phase", split.by = "updated_location")

assign("Rozowsky_seurat", integrated)
# save(Rozowsky_seurat, file = paste0(wd, "PA/PA_Data/Rozowsky.SCEAnalaysis_19.05.2022.RData"))

########## CNV Calling ##########
library(rjags)
library(infercnv)
library(Seurat)
load("PA/PA_Data/Rozowsky.SCEAnalaysis_19.05.2022.RData")
DimPlot(Rozowsky_seurat)
normal.cells.confirmed <- Cells(Rozowsky_seurat)[Rozowsky_seurat$seurat_clusters %in% c("2", "4", "6", "7", "8")]
DimPlot(Rozowsky_seurat, cells.highlight = normal.cells.confirmed)
table(Rozowsky_seurat$SingleR.Annotation, Rozowsky_seurat$seurat_clusters)

counts_matrix = GetAssayData(Rozowsky_seurat, slot="counts")
annotations_file = data.frame(Cells(Rozowsky_seurat), Rozowsky_seurat$seurat_clusters)
infercnv_obj = CreateInfercnvObject(raw_counts_matrix = counts_matrix, 
                                    annotations_file = ,
                                    # delim="\t",
                                    gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
                                    ref_group_names = c("2", "4", "6", "7", "8")) 


