# Single-cell Pilocytic Astrocytoma analysis (Reitman et al, 2019)
# Last updated: 20/12/2021

########## Run Seurat analysis on Reitman dataset ##########
library(Seurat)
library(cowplot)
scpa_counts <- read.table("/Users/jrozowsky/Documents/PMC/Data/PilocyticAstrocytoma_2019/expression/rawdata.txt", header = TRUE, row.names = 1)
scpa <- CreateSeuratObject(counts = scpa_counts, project = "PilocyticAstrocytoma", assay = "RNA", min.cells = 3, min.features = 200)
scpa
scpa@meta.data$Condition <- "PA"

scpa <- NormalizeData(scpa, normalization.method = "LogNormalize", scale.factor = 10000)
scpa <- FindVariableFeatures(scpa, selection.method = "vst", nfeatures = 4000)
scpa <- ScaleData(scpa)
scpa <- RunPCA(scpa)

VizDimLoadings(scpa, dims = 1:2, reduction = "pca")
DimPlot(scpa, reduction = "pca")
DimHeatmap(scpa, dims = 1:10, cells = 500, balanced = TRUE)

ElbowPlot(scpa, ndims = 30) # continue with 15 dims

### Find clusters ###
scpa <- FindNeighbors(scpa, reduction = "pca", dims = 1:15)

for (res in c(0.1, 0.25)) {
  scpa <- FindClusters(scpa, graph.name = "RNA_snn", resolution = res, algorithm = 1)
}
scpa <- RunTSNE(scpa, dims = 1:20)

plot_grid(ncol = 2, 
          DimPlot(scpa, reduction = "tsne", group.by = "RNA_snn_res.0.1"),
          DimPlot(scpa, reduction = "tsne", group.by = "RNA_snn_res.0.25"))

set.clust = "RNA_snn_res.0.1"
scpa <- SetIdent(scpa, value = set.clust)
DimPlot(scpa, reduction = "tsne", label = TRUE, repel = TRUE)

genes_plot <- c("APOD", "IL32", "PDGFRA", "CCL3", "SOD2", "TMEM119", "P2RY12", "S100A9", "CX3CR1")
FeaturePlot(scpa, features = genes_plot, reduction = "tsne", coord.fixed = TRUE, ncol = 3, combine = TRUE)

### Annotate clusters ###
scpa <- RenameIdents(scpa,
                     `0` = "Tumor",
                     `1` = "Microglia",
                     `2` = "Tumor",
                     `3` = "Tumor",
                     `4` = "T-cell",
                     `5` = "Microglia",
                     `6` = "Peripheral macrophage")
DimPlot(scpa, label = TRUE, repel = TRUE, label.box = TRUE, label.size = 3, label.color = "white")

marker_genes <- FindAllMarkers(scpa, logfc.threshold = 0.25, only.pos = TRUE, assay = "RNA")
top10 <- marker_genes %>% group_by(cluster) %>% top_n(-10, wt = p_val_adj)

par(mfrow = c(2, 3))
for (i in unique(top10$cluster)) {
  barplot(sort(setNames(top10$avg_log2FC, top10$gene)[top10$cluster == i], F), horiz = T,
          las = 1, main = paste0(i, " vs. rest"), border = "white", yaxs = "i")
}

### Figures ###
DotPlot(scpa, features = rev(as.character(unique(top10$gene))), group.by = set.clust, assay = "RNA") + 
  theme(axis.text.y = element_text(size = 5)) +
  coord_flip()

scpa.updated <- ScaleData(scpa, features = as.character(unique(top10$gene), assay = "RNA"))
scpa.updated$CellCluster <- scpa.updated@active.ident
DoHeatmap(scpa.updated, features = as.character(unique(top10$gene)), size = 2, draw.lines = TRUE) +
  scale_fill_continuous() + 
  theme(axis.text.y = element_text(size = 5))

########## Make Single-Cell Expression Set for MuSiC ##########
library(Biobase)
sc_cluster <- data.frame(scpa@active.ident, scpa@meta.data$orig.ident)
names(sc_cluster) <- c("Cluster", "Tumor")
sc_cluster$Tumor <- sample(sc_cluster$Tumor)

phenoData <- AnnotatedDataFrame(sc_cluster)

assayData <- as.matrix(scpa@assays$RNA@data)

# Note: have to make sure that the order of the cells in the assay and pheno data match
# check: colnames(assayData) == rownames(phenoData@data) or 
# identical(colnames(assayData), rownames(phenoData@data))
SinglePA_eset <- ExpressionSet(assayData = assayData, phenoData = phenoData)
saveRDS(SinglePA_eset, "/Users/jrozowsky/Documents/PMC/PA/PA_DataSets/PA_singleCellReitman.RData")
