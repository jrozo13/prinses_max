# Single-cell Pilocytic Astrocytoma analysis (Vladoiu et al, 2019)
# Last updated: 14/03/2022

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
DimPlot(pa.combined, reduction = "pca", group.by = "Phase")
# cancer cells should be in G1 --> don't regress out

# Determine number of dimensions to continue analysis with
# Do this by: JackStraw and Elbow plots
pa.combined <- JackStraw(pa.combined, num.replicate = 100, dims = 50)
pa.combined <- ScoreJackStraw(pa.combined, dims = 1:50)
JackStrawPlot(pa.combined, dims = 1:50)
ElbowPlot(pa.combined, ndims = 50) # continue with 30 dims

pa.combined <- RunUMAP(pa.combined, dims = 1:30, n.neighbors = 50, min.dist = 0.1)
pa.combined <- FindNeighbors(pa.combined, dims = 1:30, annoy.metric = "euclidean")
pa.combined <- FindClusters(pa.combined, resolution = 0.1, group.singletons = TRUE)
print(DimPlot(pa.combined, reduction = "umap"))

##### Functional annotation #####
library(ComplexHeatmap)
load("PA/PA_Data/SeuratObject_16.03.2022.RData")
load("PA/PA_Data/scAnalysis_16.02.2022.RData")
DimPlot(pa.combined, reduction = "umap", split.by = "orig.ident") 

uMap_cluster <- DimPlot(pa.combined, reduction = "umap", group.by = "seurat_clusters", pt.size = 1) +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  labs(title = "Pilocytic Astrocytoma",
       subtitle = "Seurat Clusters") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        plot.title = element_text(),
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))

pdf(paste0(fwd, "scUMAP.Taylor_cluster.pdf"), width = 7, height = 6)
print(uMap_cluster)
dev.off()

uMap_patient <- DimPlot(pa.combined, reduction = "umap", group.by = "orig.ident", pt.size = 1) +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  labs(title = "Pilocytic Astrocytoma",
       subtitle = "Patient") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        plot.title = element_text(),
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))

pdf(paste0(fwd, "scUMAP.Taylor_patient.pdf"), width = 7, height = 6)
print(uMap_patient)
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

pdf(paste0(fwd, "scClusterDEGs_2.pdf"), width = 10, height = 12)
draw(ht_list, heatmap_legend_side = "bottom")
dev.off()

save(markers, pa.combined_heatmap, file = paste0(wd, "PA/PA_Data/scAnalysis_16.02.2022.RData"))
save(pa.combined, file = paste0(wd, "PA/PA_Data/SeuratObject_22.03.2022.RData"))

########## SingleR ##########
library(SingleR)
library(celldex)
ref.data <- celldex::HumanPrimaryCellAtlasData(ensembl = FALSE)
# ensembl == TRUE: will use the ENSEMBL IDs
cell_predictions <- SingleR(as.SingleCellExperiment(pa.combined), 
                       ref = ref.data,
                       labels = ref.data$label.main)

cluster_predictions <- SingleR(as.SingleCellExperiment(pa.combined), 
                       ref = ref.data,
                       labels = ref.data$label.main,
                       clusters = pa.combined$seurat_clusters)

cluster_idents <- data.frame(cluster_predictions@rownames, cluster_predictions@listData$labels)
pa.combined <- RenameIdents(object = pa.combined,
                            `0` = "Macrophage",
                            `1` = "Tumor",
                            `2` = "Monocyte",
                            `3` = "Tumor",
                            `4` = "Leukeocyte",
                            `5` = "Macrophage",
                            `6` = "Macrophage")
pa.combined <- AddMetaData(pa.combined, metadata = cell_predictions@listData$labels, col.name = "SingleR_CellPred")
uMap_annot <- DimPlot(pa.combined, reduction = "umap", pt.size = 1) +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  labs(title = "Pilocytic Astrocytoma",
       subtitle = "Cluster annotations") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))
pdf(paste0(fwd, "scUMAP.Taylor_annot.pdf"), width = 7.5, height = 6)
print(uMap_annot)
dev.off()

DimPlot(pa.combined, group.by = "SingleR_CellPred")

heatmap_singleR <- plotScoreHeatmap(cell_predictions,
                 cluster_cols = FALSE,
                 annotation_col = data.frame(pa.combined@active.ident))
pdf(paste0(fwd, "heatmap_singleR.pdf"), width = 12, height = 8)
print(heatmap_singleR)
dev.off()

macrophage_cells <- subset(pa.combined, idents = "Macrophage") %>% Cells()
monocyte_cells <- subset(pa.combined, idents = "Monocyte") %>% Cells()
tam_degs <- FindMarkers(pa.combined, 
                        ident.1 = macrophage_cells, 
                        ident.2 = monocyte_cells,
                        min.pct = 0.25, max.cells.per.ident = 1000) %>% 
  data.frame %>%
  rownames_to_column(var = "Gene") %>%
  mutate(Meta = ifelse(avg_log2FC < 0 & p_val_adj < 0.001, "Enriched: Macrophage",
                       ifelse(avg_log2FC > 0 & p_val_adj < 0.001, "Enriched: Monocyte",
                              "Not significant")))
ggplot(data = tam_degs, aes(x=avg_log2FC, y=-log(p_val_adj), col = Meta)) + 
  geom_point() + 
  scale_color_manual(values=c("Red", "Blue", "Black")) +
  theme_classic() +
  ylab("-log 10 adjusted p-value") +
  xlab("log2FC") +
  theme(
    legend.title = element_blank()
  )

########## Analaysis of tumor stem cells ########## 
tumor_stem_cells <- subset(pa.combined, subset = SingleR_CellPred == c("Tissue_stem_cells")) %>% Cells()
tumor_cells <- subset(pa.combined, idents = "Tumor") %>% Cells()
setdiff(tumor_stem_cells, tumor_cells) # all TSM are tumor cells
nonTumor_stem_cells <- setdiff(tumor_cells, tumor_stem_cells)
TSM_degs <- FindMarkers(pa.combined,
            ident.1 = tumor_stem_cells, 
            ident.2 = nonTumor_stem_cells,
            #max.cells.per.ident = 200,
            min.pct = 0.25) %>% 
  data.frame %>%
  rownames_to_column(var = "Gene") %>%
  mutate(Meta = ifelse(avg_log2FC < 0 & p_val_adj < 0.001, "Down",
                       ifelse(avg_log2FC > 0 & p_val_adj < 0.001, "Up",
                              "Not significant")))
tsm_up_genes <- TSM_degs %>% 
  filter(Meta == "Up") %>%
  arrange(-avg_log2FC) %>%
  pull(Gene)
tsm_down_genes <- TSM_degs %>% 
  filter(Meta == "Down") %>%
  arrange(avg_log2FC) %>%
  pull(Gene)

tsm_volcano <- ggplot(data=TSM_degs, aes(x=avg_log2FC, y=-log(p_val_adj), col = Meta)) + 
  geom_point() + 
  scale_color_manual(values=c("Blue", "Black", "Red")) +
  theme_classic() +
  ylab("-log 10 adjusted p-value") +
  xlab("log2FC") +
  theme(panel.background = element_rect(colour = "black", size = 1),
          legend.title = element_blank(),
          axis.line = element_blank(),
          axis.title = element_text(size = 10)) +
  geom_text_repel(data = TSM_degs[TSM_degs$Gene %in% c(tsm_up_genes[1:10], tsm_down_genes[1:10]),], 
                  aes(label = Gene),
                  show.legend = FALSE
                    )
pdf(paste0(fwd, "volcano_tumorstemcell.pdf"), width = 8, height = 6)
print(tsm_volcano)
dev.off()

library(clusterProfiler)
library(org.Hs.eg.db)
library(gridExtra)
tsm_up_entrez <- AnnotationDbi::select(org.Hs.eg.db,
                    keys = tsm_up_genes,
                    columns = c("ENTREZID", "SYMBOL"),
                    keytype = "SYMBOL") %>%
  na.omit() %>%
  pull(ENTREZID)
tsm_down_entrez <- AnnotationDbi::select(org.Hs.eg.db,
                                       keys = tsm_down_genes,
                                       columns = c("ENTREZID", "SYMBOL"),
                                       keytype = "SYMBOL") %>%
  na.omit() %>%
  pull(ENTREZID)

go_tsm <- enrichGO(gene = tsm_up_entrez, 
         OrgDb = org.Hs.eg.db, 
         keyType = 'ENTREZID',
         ont = "BP",
         pAdjustMethod = "BH", 
         qvalueCutoff = 0.05, 
         readable = TRUE) %>%
  enrichplot::pairwise_termsim() %>%
  dotplot(x = "GeneRatio",
          showCategory = 20,
          font.size = 8,
          title = "Pathways enriched in tumor stem cells")
go_nonTsm <- enrichGO(gene = tsm_down_entrez, 
         OrgDb = org.Hs.eg.db, 
         keyType = 'ENTREZID',
         ont = "BP", 
         pAdjustMethod = "BH", 
         qvalueCutoff = 0.05, 
         readable = TRUE) %>%
  enrichplot::pairwise_termsim() %>%
  dotplot(x = "GeneRatio",
          showCategory = 20,
          font.size = 8,
          title = "Pathways enriched in tumor cells")

pdf(paste0(fwd, "go_tumorStemCells.pdf"), width = 12, height = 6)
grid.arrange(go_tsm, go_nonTsm, ncol = 2)
dev.off()

uMap_tsc <-  DimPlot(pa.combined, 
        cells.highlight = tumor_stem_cells, 
        cols.highlight = "red",
        sizes.highlight = 0.8,
        pt.size = 1) +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  labs(title = "Pilocytic Astrocytoma",
       subtitle = "Tumor Stem Cells") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))
  
uMap_nonTsc <- DimPlot(pa.combined,
        cells.highlight = nonTumor_stem_cells,
        cols.highlight = "blue",
        sizes.highlight = 0.8,
        pt.size = 1) +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  labs(title = "Pilocytic Astrocytoma",
       subtitle = "Other Tumor Cells") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))

pdf(paste0(fwd, "scUMAP.Taylor_tsc.pdf"), width = 12, height = 6)
grid.arrange(uMap_tsc, uMap_nonTsc, ncol = 2)
dev.off()

########## Annotation of immune cells ##########
load("PA/PA_Data/Pombo.SeuratObj_30.03.2022.RData")
assign("ref.pombo", seuratObj); rm(seuratObj)
load("PA/PA_Data/SeuratObject_22.03.2022.RData")
library(SingleR)
library(SingleCellExperiment)
library(Seurat)

immune.obj <- subset(x = pa.combined, idents = c("Leukeocyte", "Monocyte", "Macrophage"))
DefaultAssay(immune.obj) <- "integrated"

#immune.obj <- NormalizeData(immune.obj, verbose = F)
#immune.obj <- FindVariableFeatures(immune.obj, verbose = F)
immune.obj <- ScaleData(immune.obj, verbose = F)
immune.obj <- RunPCA(immune.obj, features = VariableFeatures(immune.obj))
dims.use <- 30
immune.obj <- RunUMAP(immune.obj, dims = 1:dims.use, min.dist = 0.5, verbose = F)
DimPlot(immune.obj, group.by = "orig.ident", reduction = "umapHarmony")
immune.obj <- FindNeighbors(immune.obj)
immune.obj <- FindClusters(immune.obj, group.singletons = FALSE, resolution = 0.5)
DimPlot(immune.obj, reduction = "umapHarmony")

library(harmony)
theta.use <- 1
immune.obj <- RunHarmony(immune.obj, 
                        group.by.vars = "orig.ident",
                        theta = theta.use,  
                        plot_convergence = FALSE)
immune.obj <- RunUMAP(immune.obj,
                     reduction = "harmony", 
                     dims = 1:dims.use, 
                     reduction.name = "umapHarmony",
                     reduction.key = "umapHarmony")

DimPlot(immune.obj,
        group.by = "orig.ident",
        reduction = "umapHarmony") +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        plot.title = element_text(),
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))

DimPlot(immune.obj,
        reduction = "umapHarmony",
        group.by = "SingleR_CellPred") +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        plot.title = element_text(),
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))

ref.pombo <- as.SingleCellExperiment(ref.pombo)
immune.obj_sce <- as.SingleCellExperiment(immune.obj)

immune.cell_predictions <- SingleR(test = immune.obj_sce,
                                   ref = ref.pombo,
                                   labels = ref.pombo$cluster)

immune.cluster_predictions <- SingleR(test = immune.obj_sce,
                                   ref = ref.pombo,
                                   clusters = immune.obj_sce$seurat_clusters,
                                   labels = ref.pombo$cluster)

save(immune.cell_predictions, immune.cluster_predictions, file = paste0(wd, "PA/PA_Data/SingleRAnnotation_Taylor.RData"))

immune.obj <- AddMetaData(immune.obj, immune.cell_predictions$labels, col.name = "SingleR.Pombo_Cell")


singleR_clusterpred <- data.frame(SeuratCluster = immune.obj$seurat_clusters %>% levels(),
                                  PomboCellType = immune.cluster_predictions$labels)

immuneCell_clusterPred <- data.frame(SeuratCluster = immune.obj$seurat_clusters) %>% 
  rownames_to_column(var = "Cell")
immuneCell_clusterPred <- merge(x = immuneCell_clusterPred, 
                                y = singleR_clusterpred, 
                                by.x = "SeuratCluster", 
                                by.y = "SeuratCluster") 
immuneCell_clusterPred <- immuneCell_clusterPred %>% arrange(Cell)
(immuneCell_clusterPred$SeuratCluster == immune.obj$seurat_clusters) %>% table()
immune.obj <- AddMetaData(immune.obj, immuneCell_clusterPred$PomboCellType, col.name = "SingleR.Pombo_Cluster")

uMapPA_pomboAnnot <- DimPlot(immune.obj,
        reduction = "umapHarmony",
        group.by = "SingleR.Pombo_Cluster",
        pt.size = 1) +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  labs(title = "Pilocytic Astrocytoma: immune microenvironment",
       subtitle = "SingleR annotation: Pombo-Antunes, 2021") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        plot.title = element_text(),
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))

pdf(paste0(fwd, "scUMAP.TaylorPA_pomboAnnot.pdf"), width = 7.5, height = 6)
print(uMapPA_pomboAnnot)
dev.off()

save(immune.obj, file = paste0(wd, "PA/PA_Data/TaylorImmune.SeuratObj_06.04.2022.RData"))

iummne.markers <- FindAllMarkers(immune.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5_immune <- leuk.markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
FeaturePlot(immune.obj, features = "S100A6")
DimPlot(immune.obj)
scDEG_heatmap <- DoHeatmap(immune.obj, features = top5_immune$gene) + NoLegend()

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

