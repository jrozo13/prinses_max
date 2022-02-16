# Single-cell Pilocytic Astrocytoma analysis (Vladoiu et al, 2019)

wd <- "/Users/jrozowsky/Library/Mobile Documents/com~apple~CloudDocs/Documents/PMC/"
setwd(wd)
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
        features = c("nCount_RNA", "nFeature_RNA", "percent.mito"), nCol = 3,  point.size.use = 0.1)
VlnPlot(pa.combined, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"), ncol = 3)


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
