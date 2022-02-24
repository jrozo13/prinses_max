# Single-cell Pilocytic Astrocytoma analysis (Vladoiu et al, 2019)
# Last updated: 17/02/2022

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

########## Load data ##########
library(Seurat)
# Compile all patient IDs. Make Seurat objects with for-loop
patient_ids <- list.files(path = paste0(wd, "Data/Vladoiu_2019/"))
for (patient in patient_ids) {
  file_dir <- paste0(wd, "Data/Vladoiu_2019/", patient, "/")
  counts_name <- paste0(patient, "_counts")
  print(file_dir)
  print(obj_name)
  
  counts <- Read10X(data.dir = file_dir)
  assign(counts_name, counts)
    
  object <- CreateSeuratObject(counts = counts, project = patient, assay = "RNA")
  assign(obj_name, object)
}

### Analyze sequencing metrics of single cells
counts_per_cell <- Matrix::colSums(PA0706_counts)
counts_per_gene <- Matrix::rowSums(PA0706_counts)
genes_per_cell <- Matrix::colSums(PA0706_counts>0)
hist(log10(counts_per_cell+1),main='counts per cell',col='wheat')
hist(log10(genes_per_cell+1), main='genes per cell', col='wheat')
plot(counts_per_cell, genes_per_cell, log='xy', col='wheat')
plot(sort(genes_per_cell), xlab='cell', log='y', main='genes per cell (ordered)')

pa.combined <- merge(x = PA0706_ojb, y = PA2406_ojb,
                     add.cell.ids = c(patient_ids[1], patient_ids[2]), 
                     project = "PiloAstro")

########## Pre-process object ##########
mito.genes <- grep(pattern = "^MT-", x = rownames(x = pa.combined@assays$RNA@counts), value = TRUE)
percent.mito <- Matrix::colSums(pa.combined@assays$RNA@counts[mito.genes, ])/Matrix::colSums(pa.combined@assays$RNA@counts)
pa.combined <- AddMetaData(object = pa.combined, metadata = percent.mito, col.name = "percent.mito")
hist(pa.combined$percent.mito)
table(pa.combined$percent.mito < 0.05) # 4939 cells have >5% mitochondrial gene expression
head(pa.combined@meta.data)
VlnPlot(object = pa.combined,
        features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), 
        ncol = 3, group.by = "orig.ident", pt.size = 0)

# Filter cells that have < 200 and > 5000 genes expressed, and > 10% mitochondria reads
pa.combined <- subset(x = pa.combined, 
                      subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mito < 10)
# we are left with 15823 cells

##### Split combined data-set by patient #####
# Normalize and find variable features separately
pa.separated <- SplitObject(pa.combined, split.by = "orig.ident")
pa.separated <- lapply(X = pa.separated, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 1e4)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Find genes to integrate samples by

keep_features <- c("TMEM119", "P2RY12", "CD163", "CCL2", "ITGAL", "GFAP", "CXCL11", "MOG", "PLP1", "TOX", "PDCD1", "LAG3", "TIGIT")

features <- SelectIntegrationFeatures(object.list = pa.separated)
pa.anchors <- FindIntegrationAnchors(object.list = pa.separated, anchor.features = c(features, keep_features))
pa.combined <- IntegrateData(anchorset = pa.anchors)

##### Process Integrated Data Set #####
DefaultAssay(pa.combined) <- "integrated"
pa.combined <- FindVariableFeatures(object = pa.combined,
                                    selection.method = "vst",
                                    nfeatures = 2000)
pa.combined <- ScaleData(pa.combined, do.center = TRUE)

########## Clustering ##########
pa.combined <- RunPCA(pa.combined, features = VariableFeatures(object = pa.combined))
VizDimLoadings(pa.combined, dims = 1:2, reduction = "pca")
DimPlot(pa.combined, reduction = "pca", group.by = "orig.ident")

# Determine number of dimensions to continue analysis with
# Do this by: JackStraw and Elbow plots
pa.combined <- JackStraw(pa.combined, num.replicate = 100, dims = 50)
pa.combined <- ScoreJackStraw(pa.combined, dims = 1:50)
JackStrawPlot(pa.combined, dims = 1:50)
ElbowPlot(pa.combined, ndims = 50) # continue with 30 dims

pa.combined <- RunUMAP(pa.combined, dims = 1:30, n.neighbors = 30, min.dist = 0.5)
pa.combined <- FindNeighbors(pa.combined, dims = 1:30, annoy.metric = "euclidean")
pa.combined <- FindClusters(pa.combined, resolution = 0.3, group.singletons = TRUE)
print(DimPlot(pa.combined, reduction = "umap"))
save(pa.combined, file = paste0(wd, "PA/PA_Data/SeuratObject_16.02.2022.RData"))

##### Functional annotation #####
library(ComplexHeatmap)
load("PA/PA_Data/SeuratObject_16.02.2022.RData")
DimPlot(pa.combined, reduction = "umap", split.by = "orig.ident") 

pdf(paste0(fwd, "scUMAP_cluster.pdf"), width = 6, height = 6)
print(DimPlot(pa.combined, reduction = "umap"))
dev.off()

pdf(paste0(fwd, "scUMAP_patient.pdf"), width = 6, height = 6)
print(DimPlot(pa.combined, reduction = "umap", group.by = "orig.ident"))
dev.off()

markers <- FindAllMarkers(pa.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
scDEG_heatmap <- DoHeatmap(pa.combined, features = top5$gene) + NoLegend()

pa.combined_heatmap <- AverageExpression(pa.combined, 
                       assays = "integrated",
                       features = top5$gene,
                       group.by = "seurat_clusters",
                       return.seurat = TRUE)

col_fun = colorRamp2(c(-3, 0, 3), c("navy", "white", "red"))
ht_list <- Heatmap(pa.combined_heatmap@assays$integrated@scale.data,
        cluster_rows = FALSE,
        cluster_columns = FALSE, 
        border = "black",
        column_split = c(0:10),
        col = col_fun,
        show_column_names = FALSE, 
        heatmap_legend_param = list(
          title = "Relative expression",
          direction = "horizontal",
          at = c(-3, 0, 3),
          border = "black",
          legend_width = unit(6, "cm"),
          title_position = "topcenter"
          ))

pdf(paste0(fwd, "scClusterDEGs.pdf"), width = 6, height = 12)
draw(ht_list, heatmap_legend_side = "bottom")
dev.off()

save(markers, pa.combined_heatmap, file = paste0(wd, "PA/PA_Data/scAnalysis_16.02.2022.RData"))

##### SingleR #####
Seurat::FindSubCluster()

##### SingleR #####
library(SingleR)
library(celldex)


########## CNV ##########
# library(infercnv)
# counts_matrix = as.matrix(pa.combined@assays$RNA@counts[,colnames(pa.combined)])


########## Make Single-Cell Expression Set for MuSiC ##########
library(Biobase)
sc_cluster <- data.frame(pa.combined@active.ident, pa.combined@meta.data$orig.ident)
names(sc_cluster) <- c("Cluster", "Sample")
table(sc_cluster$Cluster, sc_cluster$Sample)

phenoData <- AnnotatedDataFrame(sc_cluster)

assayData <- as.matrix(pa.combined@assays$RNA@data)

# Note: have to make sure that the order of the cells in the assay and pheno data match
# check: colnames(assayData) == rownames(phenoData@data) or 
# identical(colnames(assayData), rownames(phenoData@data))
Vladoiu_eset <- ExpressionSet(assayData = assayData, phenoData = phenoData)
save(Vladoiu_eset,
     file = paste0(wd, "PA/Analysis/Deconvolution/Vladoiu_scEset.RData"))

########## Make Single-Cell Expression Set for MuSiC ##########
load("Data/PMC/BulkPA_eset.RData")
