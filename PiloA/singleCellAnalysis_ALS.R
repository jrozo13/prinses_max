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

########## Load data ##########
# compile all patient IDs
patient_info <- read.delim(file = "Data/ALS_2022/SCPCP000002/single_cell_metadata.tsv") %>%
  filter(diagnosis == "Pilocytic astrocytoma")
patient_ids <- patient_info %>%
  pull(scpca_sample_id)
library_ids <- patient_info %>%
  pull(scpca_library_id)
patient_dict <- data.frame(patient_ids, library_ids)

st_loc <- c("Basal ganglia", "Left mesial temporal", "Supracellar", "Thalamic", "Thalamus")
pf_loc <- c("Cerebellar", "Posterior fossa")
patient_info <- patient_info %>%
  mutate(updated_location = ifelse(tissue_location %in% st_loc, "Supratentorial",
                                   ifelse(tissue_location %in% pf_loc, "Posterior fossa", "Spinal")))

# make Seurat objects for each sample
library(Seurat)
library(SingleCellExperiment)
for (patient in patient_ids) {
  library <-  patient_dict$library_ids[patient_dict$patient_ids==patient]
  object <- readRDS(paste0(wd, "Data/ALS_2022/SCPCP000002/", patient, "/", library, "_filtered.rds"))
  
  # convert from ensmbl ID to gene symbol
  rownames(object) <- object@rowRanges@elementMetadata$gene_symbol
  
  # remove rows where ensmbl ID didnt map to gene symbol
  object <- object[!is.na(rownames(object)),]
  
  seuratObj <- CreateSeuratObject(counts = assay(object), project = patient, assay = "RNA")
  
  print(paste0(patient, "_sce"))
  assign(paste0(patient, "_sce"), seuratObj)
  rm(seuratObj)
  rm(object)
}

# merge all individual Seurat objects
combined <- merge(get(paste0(patient_ids[1], "_sce")),
                  y = mget(paste0(patient_ids[-1], "_sce")),
                  add.cell.ids = c(patient_ids),
                  project = "ALS_PiloAstro")
rm(list = paste0(patient_ids, "_sce"))

pf.sample <- patient_info$scpca_sample_id[patient_info$updated_location == "Posterior fossa"]
st.sample <- patient_info$scpca_sample_id[patient_info$updated_location == "Supratentorial"]
spine.sample <- patient_info$scpca_sample_id[patient_info$updated_location == "Spinal"]

meta.data.keep <- cbind(combined@meta.data)
meta.data.keep <- meta.data.keep %>% 
  as.data.frame() %>%
  mutate(updated_location = ifelse(meta.data.keep$orig.ident %in% pf.sample, "Posterior fossa",
                                   ifelse(meta.data.keep$orig.ident %in% st.sample, "Supratentorial", 
                                          ifelse(meta.data.keep$orig.ident %in% spine.sample, "Spinal", ""))))

combined <- AddMetaData(combined, meta.data.keep)

########## Quality Control ########## 
library(scater)
# change this to the object of interest
seurat.object = combined

sce <- as.SingleCellExperiment(x = seurat.object)
sce@colData$nCount_RNA <- NULL
sce@colData$nFeature_RNA <- NULL

# some quality results are in the seurat object metadata, but want to calculate on our own
mito.genes <- grep("^MT-", rownames(sce), value = TRUE) # 37 mitochondrial genes
sce$percent.mito <- (Matrix::colSums(counts(sce)[mito.genes, ])*100)/Matrix::colSums(counts(sce)) 
sce$nGene <- apply(counts(sce), 2, function(x) length(x[x > 0])) # number of expressed genes
sce$nUMI <- apply(counts(sce), 2, sum) # total UMI counts (library size)
sce$staticNr <- 1
dim(colData(sce))
colnames(colData(sce))

summary(sce$nUMI)
ggplot(mapping = aes(x = sce$nUMI)) + 
  geom_density(fill = "lightblue") +
  geom_vline(xintercept = 1000, color = "red") +
  labs(x = "Counts per cell") +
  xlim(0, 30000)
sum(sce$nUMI < 1000)

summary(sce$nGene)
ggplot(mapping = aes(x = sce$nGene)) + 
  geom_density(fill = "lightblue") +
  geom_vline(xintercept = 500, color = "red") +
  labs(x = "Number of genes expressed")
sum(sce$nGene < 500)

ggplot(mapping = aes(x = sce$percent.mito)) +
  geom_density(fill = "lightblue") + 
  labs(x = "Mitchondrial fraction") +
  geom_vline(xintercept = 20, color = "red")
ggplot(mapping = aes(x = sce$percent.mito)) +
  geom_density(fill = "lightblue") + 
  labs(x = "Mitchondrial fraction") +
  geom_vline(xintercept = 25, color = "red")
ggplot(mapping = aes(x = sce$percent.mito)) +
  geom_density(fill = "lightblue") + 
  labs(x = "Mitchondrial fraction") +
  geom_vline(xintercept = 30, color = "red")

# Set Thresholds for qualtiy control
genes_exp.check = 300
total_counts.check = 1000
mito_fraction.check = 25

qc_df <- data.frame(barcode = Cells(sce),
                    genes_exp = sce$nGene,
                    total_counts = sce$nUMI,
                    mito_fraction = sce$percent.mito)
qc_df <- qc_df %>% 
  mutate(filtered = 
           ifelse(genes_exp >= genes_exp.check & 
                    total_counts >= total_counts.check & 
                    mito_fraction <= mito_fraction.check, "Passed", "Failed"))

ggplot(qc_df, aes (x = total_counts,
                   y = genes_exp, 
                   color = filtered)) +
  geom_point(alpha = 0.5) +
  geom_vline(xintercept = total_counts.check, color = "red") +
  geom_hline(yintercept = genes_exp.check, color = "red") +
  labs(x = "Total Count",
       y = "Number of Genes Expressed",
       color = paste0("Mitochondrial\nFraction: ", mito_fraction.check, "%")) + 
  theme_bw()
ggplot(qc_df, aes (x = total_counts,
                   y = genes_exp, 
                   color = mito_fraction)) +
  geom_point(alpha = 0.5) +
  scale_color_viridis_c() +
  geom_vline(xintercept = total_counts.check, color = "red") +
  geom_hline(yintercept = genes_exp.check, color = "red") +
  labs(x = "Total Count",
       y = "Number of Genes Expressed",
       color = paste0("Mitochondrial\nFraction")) + 
  theme_bw()

filtered_samples <- qc_df %>%
  dplyr::filter(total_counts >= total_counts.check,
                genes_exp >= genes_exp.check,
                mito_fraction <= mito_fraction.check)

# Remove outliers (101 cells which had lowly expressed genes)
sce_clean <- sce[,rownames(filtered_samples)]
sce_clean@colData <- sce_clean@colData[Cells(sce_clean),]

# save(sce, sce_clean, file = paste0(wd, "PA/PA_Data/ALS.SingleCellExpObject_25.04.2022.RData"))

########## Make Seurat Object ##########
# load(paste0(wd, "PA/PA_Data/ALS.SingleCellExpObject_25.04.2022.RData"))
counts <- counts(sce_clean)
rownames(counts) <- rownames(sce_clean)
colnames(counts) <- colnames(sce_clean)
allcell.obj <- CreateSeuratObject(counts = counts)
allcell.obj$updated_location <- sce_clean$updated_location
allcell.obj$nUMI <- sce_clean$nUMI
allcell.obj$nGene <- sce_clean$nGene
allcell.obj$percent.mito <- sce_clean$percent.mito
allcell.obj$nCount_RNA <- NULL
allcell.obj$nFeature_RNA <- NULL

library(Seurat)
object.list <- SplitObject(allcell.obj, split.by = "orig.ident")
for (i in 1:length(object.list)) {
  print(paste0("index: ", i, "; nCells: ", length(Cells(object.list[[i]]))))
}
# remove samples with < 30 cells after filtering
object.list[4] <- NULL

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
integrated <- RunTSNE(integrated, dims = 1:dims.use, verbose = F, )
# assign(paste0("integrated_mito", mito_fraction.check), integrated)
DimPlot(integrated, reduction = "tsne", group.by = "orig.ident")
DimPlot(integrated, reduction = "tsne", group.by = "updated_location")

# We do not see separation by patient or batches
FeaturePlot(integrated,
            reduction = "tsne", 
            features = c("nUMI", "nGene", "percent.mito"))
FeaturePlot(integrated,
            reduction = "tsne", 
            features = c("MKI67"))

##### Plot difference in cells filtered out by mitochondrial read % #####
filteredCells <- setdiff(Cells(integrated_mito30), Cells(integrated_mito20))
DimPlot(integrated_mito30,
        pt.size = 0.75, 
        cells.highlight = filteredCells, 
        sizes.highlight = 0.1,
        reduction = "umap") +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  labs(title = "Pilocytic Astrocytoma",
       subtitle = "Alex's Lemonade Stand, 2022\n Filtered cells at 20% vs 30%") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))


########## Clustering ##########
integrated <- FindNeighbors(integrated,
                           dims = 1:dims.use,
                           annoy.metric = "euclidean")
integrated <- FindClusters(integrated, 
                           resolution = 1, 
                           group.singletons = FALSE,
                           algorithm = 1)
DimPlot(integrated,
        pt.size = 1,
        reduction = "tsne")

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
DimPlot(integrated, reduction = "tsne", group.by = "SingleR.Annotation", pt.size = 1.5)
DimPlot(integrated, reduction = "tsne", group.by = "integrated_snn_res.0.5")

glioma.cells <- colnames(integrated)[integrated$SingleR.Annotation == "Astrocyte"]
immune.cells <- colnames(integrated)[!integrated$SingleR.Annotation == "Astrocyte"]

immune.obj <- subset(x = integrated, cells = immune.cells)
glioma.obj <- subset(x = integrated, cells = glioma.cells)

# save(integrated,
#      immune.obj,
#      glioma.obj,
#      file = paste0(wd, "PA/PA_Data/ALS.IntegratedObject_18.05.2022.RData"))

########## Annotation of immune cells ##########
load("PA/PA_Data/Pombo.SeuratObj_30.03.2022.RData")
load("PA/PA_Data/ALS.ImmuneObject_28.04.2022.RData")
assign("ref.pombo", seuratObj); rm(seuratObj)
library(SingleR)
library(SingleCellExperiment)
library(Seurat)

DefaultAssay(immune.obj) <- "integrated"
immune.obj <- ScaleData(immune.obj, verbose = F)
immune.obj <- RunPCA(immune.obj, features = VariableFeatures(immune.obj))
dims.use <- 30
immune.obj <- RunUMAP(immune.obj, dims = 1:dims.use, min.dist = 0.5, verbose = F)
DimPlot(immune.obj, group.by = "updated_location", reduction = "umap")
DimPlot(immune.obj, group.by = "orig.ident", reduction = "umap")
# We do not see inter-patient heterogeneity in the cell types

# Under-cluster cells and annotate using SingleR and Pombo annotations
immune.obj <- FindNeighbors(immune.obj)
immune.obj <- FindClusters(immune.obj, resolution = 1)
DimPlot(immune.obj, group.by = "seurat_clusters")

### Compare our annotations to pombo and singleR annotations
ref.pombo <- as.SingleCellExperiment(ref.pombo)

# Update these parameters based on reference scSeq set
ref.set = ref.pombo
col.name = "SingleR.Pombo"
immune.cluster_predictions <- SingleR(as.SingleCellExperiment(immune.obj), 
                               ref = ref.set,
                               labels = ref.set$cluster, # udpate this line based on reference set
                               clusters = immune.obj$seurat_clusters)
immune.cluster_idents <- data.frame(immune.cluster_predictions@rownames, immune.cluster_predictions@listData$labels)
immune.cell_idents <- merge(x = immune.obj@meta.data %>% rownames_to_column() %>% select(rowname, seurat_clusters), 
                            y = immune.cluster_idents,
                            by.x = "seurat_clusters",
                            by.y = "immune.cluster_predictions.rownames") %>%
  column_to_rownames(var = "rowname")
immune.cell_idents <- immune.cell_idents[rownames(immune.obj@meta.data),]
colnames(immune.cell_idents) <- c("seurat_clusters", "SingleR.Pombo")
immune.cell_idents <- immune.cell_idents %>% mutate(CellType = ifelse(SingleR.Pombo %in% c("TAM 1", "TAM 2"), "TAM",
                           ifelse(SingleR.Pombo == "NK cells", "T cells", SingleR.Pombo)))
  
immune.obj <- AddMetaData(immune.obj, metadata = immune.cell_idents$CellType, col.name = col.name)
DimPlot(immune.obj, reduction = "umap", group.by = col.name)

scUMAP.ALS_immuneMicro <- DimPlot(immune.obj,
        reduction = "umap",
        group.by = "SingleR.Pombo",
        #split.by = "updated_location",
        pt.size = 1) +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        legend.position = "none",
        plot.title = element_blank(),
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))
pdf(paste0(fwd, "scUMAP.ALS_immuneMicro.pdf"), width = 4, height = 4)
print(scUMAP.ALS_immuneMicro)
dev.off()

Idents(immune.obj) <- "SingleR.Pombo"
markers <- FindAllMarkers(immune.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top7 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 7, wt = avg_log2FC)

immune.cluterHeatmap <- AverageExpression(immune.obj, 
                                       assays = "integrated",
                                       features = c(top7$gene),
                                       group.by = "SingleR.Pombo",
                                       return.seurat = TRUE)

col_fun = colorRamp2(c(-1.5, 0, 1.5), c("navy", "white", "red"))
Heatmap(immune.cluterHeatmap@assays$integrated@scale.data,
        cluster_rows = FALSE,
        cluster_columns = FALSE, 
        border = "black",
        column_split = c(0:5),
        col = col_fun,
        show_column_names = FALSE, 
        heatmap_legend_param = list(
          title = "Relative expression",
          direction = "horizontal",
          at = c(-2, 0, 2),
          border = "black",
          legend_width = unit(6, "cm"),
          title_position = "topcenter"
        ))

pdf(paste0(fwd, "test.pdf"), width = 3, height = 6)
draw(ht_list, heatmap_legend_side = "bottom")
dev.off()


save(immune.obj,
     immune.cluterHeatmap,
     file = paste0(wd, "PA/PA_Data/ALS.ImmuneObject_28.04.2022.RData"))

cell_annotations <- immune.obj@meta.data %>%
  select("orig.ident", "SingleR.Pombo")
cell_annotations$SingleR.Pombo <- str_replace(cell_annotations$SingleR.Pombo, pattern = "NK cells", replacement = "T/NK cells")
cell_annotations <- table(cell_annotations$SingleR.Pombo, cell_annotations$orig.ident) %>% 
  data.frame()
colnames(cell_annotations) <- c("Cell", "Sample", "Count")
ggplot(cell_annotations, aes(x = Sample, y = Count, fill = Cell)) + 
  geom_bar(stat = "identity", position="fill")

########## Annotation of glima associated macrophages ##########
library(SingleCellExperiment)
library(Seurat)
library(ComplexHeatmap)

Idents(immune.obj) <- "SingleR.Pombo"
gam.obj <- subset(x = immune.obj, idents = setdiff(immune.obj$SingleR.Pombo %>% unique(), c("NK cells", "B cells")))
# see if there is a way to select specific cells, not clusters...
DefaultAssay(gam.obj) <- "integrated"
gam.obj <- ScaleData(gam.obj, verbose = F)
gam.obj <- RunPCA(gam.obj, features = VariableFeatures(gam.obj))
dims.use <- 30
gam.obj <- RunUMAP(gam.obj, dims = 1:dims.use, min.dist = 0.5, verbose = F)
DimPlot(gam.obj, group.by = "updated_location", reduction = "umap")
DimPlot(gam.obj, group.by = "orig.ident", reduction = "umap")

# Under-cluster cells and annotate using SingleR and Pombo annotations
gam.obj <- FindNeighbors(gam.obj)
gam.obj <- FindClusters(gam.obj, group.singletons = TRUE, resolution = 0.1)
DimPlot(gam.obj, group.by = "seurat_clusters")
DimPlot(gam.obj, group.by = "updated_location")

FeaturePlot(gam.obj, features = "TMEM119")

load("PA/PA_Data/ALS.GAMObject_28.04.2022.RData")
als.umap_gams <- DimPlot(gam.obj,
                                      reduction = "umap",
                                      group.by = "seurat_clusters",
                                      pt.size = 1) +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        legend.position = "none",
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))
als.umap_gams.loc <- DimPlot(gam.obj,
                         reduction = "umap",
                         group.by = "updated_location",
                         pt.size = 1) +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        #legend.position = "none",
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))
pdf(paste0(fwd, "scUMAP.ALS_GAMs.pdf"), width = 4, height = 4)
print(als.umap_gams)
dev.off()

pdf(paste0(fwd, "scUMAP.ALS_GAMs.loc.pdf"), width = 5.5, height = 4)
print(als.umap_gams.loc)
dev.off()

dam_genes <- c("SPP1", "")

markers <- FindAllMarkers(gam.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top7 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 7, wt = avg_log2FC)

gam.cluterHeatmap <- AverageExpression(gam.obj, 
                                   assays = "integrated",
                                   features = c(top5$gene),
                                   group.by = "seurat_clusters",
                                   return.seurat = TRUE)

col_fun = colorRamp2(c(-1.5, 0, 1.5), c("navy", "white", "red"))
ht_list <- Heatmap(gam.cluterHeatmap@assays$integrated@scale.data,
        cluster_rows = FALSE,
        cluster_columns = FALSE, 
        border = "black",
        column_split = c(0:5),
        col = col_fun,
        show_column_names = FALSE, 
        heatmap_legend_param = list(
          title = "Relative expression",
          direction = "horizontal",
          at = c(-2, 0, 2),
          border = "black",
          legend_width = unit(6, "cm"),
          title_position = "topcenter"
        ))

pdf(paste0(fwd, "scCluster.ALS_GAMs.pdf"), width = 3, height = 6)
draw(ht_list, heatmap_legend_side = "bottom")
dev.off()

library(org.Hs.eg.db)
library(enrichplot)
library(clusterProfiler)
for (c in c(0:5)) {
  clusterGenes <- markers %>%
    filter(cluster == c & p_val_adj < 0.01) %>%
    arrange(-avg_log2FC) %>%
    pull(gene)
  rankGenes <- markers %>%
    filter(cluster == c & p_val_adj < 0.01) %>%
    arrange(-avg_log2FC) %>%
    pull(avg_log2FC)
  geneWithRanks <- data.frame(GeneName = clusterGenes, GeneRank = rankGenes)
  # write.table(x = geneWithRanks,
  #       file = paste0(wd, "PA/PA_Data/scGeneSets/PA_TAM.C", c, "_11.05.2022.txt"),
  #       sep = "\t",
  #       row.names = FALSE)
  
  # write(x = clusterGenes,
  #       file = paste0(wd, "PA/PA_Data/scGeneSets/PA_TAM.gProf.C", c, "_11.05.2022.txt"),
  #       sep = "\t")
  
  clusterGenes <- AnnotationDbi::select(org.Hs.eg.db,
                                        keys = clusterGenes,
                                        columns = c("ENTREZID", "SYMBOL"),
                                        keytype = "SYMBOL") %>%
    na.omit() %>%
    pull(ENTREZID)
  goPlot <- enrichGO(gene = clusterGenes,
                     OrgDb = org.Hs.eg.db,
                     keyType = 'ENTREZID',
                     ont = "BP",
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.001,
                     readable = TRUE) %>%
    enrichplot::pairwise_termsim() %>%
    dotplot(x = "GeneRatio",
            showCategory = 20,
            font.size = 8,
            title = paste0("Pathways enriched in cluster ", c))
  assign(paste0("goPlot_TAM.C", c), goPlot)
}
# remove quotation marks in .txt files
# load into GSEA
# look for enrichment of M0 v M1 phenotype
# https://www.gsea-msigdb.org/gsea/msigdb/cards/GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_UP.html
# https://www.gsea-msigdb.org/gsea/msigdb/cards/GSE5099_CLASSICAL_M1_VS_ALTERNATIVE_M2_MACROPHAGE_DN.html

Idents(ref.pombo) <- "cluster"
tam_markers <- FindMarkers(ref.pombo, ident.1 = "TAM 1", ident.2 = "TAM 2")
pombo_markers <- FindAllMarkers(ref.pombo)
pombo_markers <- pombo_markers %>%
  filter(gene %in% rownames(gam.obj))
nGene <- 15
mo_tam <- pombo_markers %>%
  filter(cluster == "TAM 1") %>%
  arrange(-avg_log2FC) %>%
  pull(gene) %>% 
  .[1:nGene]
mg_tam <- pombo_markers %>%
  filter(cluster == "TAM 2") %>%
  arrange(-avg_log2FC) %>%
  pull(gene) %>% 
  .[1:nGene]
mg_tam <- c(mg_tam, "TMEM119")
prolif_tam <- pombo_markers %>%
  filter(cluster == "prol. TAM") %>%
  arrange(-avg_log2FC) %>%
  pull(gene) %>% 
  .[1:nGene]
mono_tam <- pombo_markers %>%
  filter(cluster == "Monocytes") %>%
  arrange(-avg_log2FC) %>%
  pull(gene) %>% 
  .[1:nGene]
dc_tam <- pombo_markers %>%
  filter(cluster == "DC") %>%
  arrange(-avg_log2FC) %>%
  pull(gene) %>% 
  .[1:nGene]

# Plotting module scores
gam.obj <- AddModuleScore(gam.obj,
                       features = list(mg_tam),
                       name ="ModuleScore_mgTAM",
                       assay = "RNA")
gam.obj <- AddModuleScore(gam.obj,
                          features = list(mo_tam),
                          name ="ModuleScore_moTAM",
                          assay = "RNA")
gam.obj <- AddModuleScore(gam.obj,
                          features = list(prolif_tam),
                          name ="ModuleScore_proTAM",
                          assay = "RNA")
gam.obj <- AddModuleScore(gam.obj,
                          features = list(mono_tam),
                          name ="ModuleScore_monocyte",
                          assay = "RNA")
gam.obj <- AddModuleScore(gam.obj,
                          features = list(dc_tam),
                          name ="ModuleScore_DC",
                          assay = "RNA")

als.umap_gams.mgModuleScore <- FeaturePlot(gam.obj,
                                 features = c("ModuleScore_mgTAM1")) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu"))) +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        legend.position = "none",
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))
pdf(paste0(fwd, "scUMAP.ALS_GAMs.mgModuleScore.pdf"), width = 2, height = 2)
print(als.umap_gams.mgModuleScore)
dev.off()

als.umap_gams.moModuleScore <- FeaturePlot(gam.obj,
                                                  features = c("ModuleScore_moTAM1")) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu"))) +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        legend.position = "none",
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))
pdf(paste0(fwd, "scUMAP.ALS_GAMs.moTamModuleScore.pdf"), width = 2, height = 2)
print(als.umap_gams.moModuleScore)
dev.off()

als.umap_gams.prolifTamModuleScore <- FeaturePlot(gam.obj,
                                           features = c("ModuleScore_proTAM1")) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu"))) +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        legend.position = "none",
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))
pdf(paste0(fwd, "scUMAP.ALS_GAMs.prolifTamModuleScore.pdf"), width = 2, height = 2)
print(als.umap_gams.prolifTamModuleScore)
dev.off()



FeaturePlot(gam.obj,
            features = c("ModuleScore_proTAM1")) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu")))

save(gam.obj, gam.cluterHeatmap, markers, pombo_markers, file = paste0(wd, "PA/PA_Data/ALS.GAMObject_28.04.2022.RData"))

########## Annotation of lymphocytes ##########
library(SingleCellExperiment)
library(Seurat)
library(ComplexHeatmap)

Idents(immune.obj) <- "SingleR.Pombo"
immune.obj$SingleR.Annotation %>% table()
immune.obj$SingleR.Pombo %>% table()
# the number of T-cell + NK-cell is same across both annotations --> check
# the number of B-cell is the same across both annotations --> check

lym.obj <- subset(x = immune.obj, idents = c("NK cells"))
# see if there is a way to select specific cells, not clusters...
DefaultAssay(lym.obj) <- "integrated"
lym.obj <- ScaleData(lym.obj, verbose = F)
lym.obj <- RunPCA(lym.obj, features = VariableFeatures(lym.obj))
dims.use <- 30
lym.obj <- RunUMAP(lym.obj, dims = 1:dims.use, min.dist = 0.5, verbose = F)
DimPlot(lym.obj, group.by = "updated_location", reduction = "umap")
DimPlot(lym.obj, group.by = "orig.ident", reduction = "umap")
DimPlot(lym.obj, group.by = "SingleR.Annotation", reduction = "umap")

myeloid_genes <- c("C1QA", "CD163", "C1QB", "C1QC")
t_genes <- c("CD3D", "CD3E", "CD3G", "CD247")

lym.obj <- AddModuleScore(lym.obj,
                          features = list(myeloid_genes),
                          name = "ModuleScore_myeloid",
                          assay = "RNA")
FeaturePlot(lym.obj, features = "ModuleScore_myeloid1")

lym.obj <- AddModuleScore(lym.obj,
                          features = list(t_genes),
                          name = "ModuleScore_TCell",
                          assay = "RNA")
FeaturePlot(lym.obj, features = "ModuleScore_TCell1")

remove.myeloid <- colnames(lym.obj)[lym.obj$ModuleScore_myeloid1 > 0]
remove.both <- colnames(lym.obj)[lym.obj$ModuleScore_TCell1 > 0.5 & lym.obj$ModuleScore_myeloid1 > 0.5]

lym.obj <- lym.obj[,!colnames(lym.obj) %in% c(remove.myeloid, remove.both)]
DimPlot(lym.obj, group.by = "SingleR.Annotation", reduction = "umap")

ggplot(data = lym.obj@meta.data, aes(x = ModuleScore_TCell1, y = ModuleScore_myeloid1)) + geom_point()
# Cluster 0 expresses high levels of myeloid and t cell related genes --> remove cluster 0

### Proceed with filtered cells
DefaultAssay(lym.obj) <- "integrated"
lym.obj <- ScaleData(lym.obj, verbose = F)
lym.obj <- RunPCA(lym.obj, features = VariableFeatures(lym.obj))
dims.use <- 30
lym.obj <- RunUMAP(lym.obj, dims = 1:dims.use, min.dist = 1, verbose = F)
lym.obj <- FindNeighbors(lym.obj)
lym.obj <- FindClusters(lym.obj, group.singletons = TRUE, resolution = 0.1)
DimPlot(lym.obj, group.by = "seurat_clusters", reduction = "umap")

nk_genes <- c("KIR2DL1", "KIR2DL3", "KIR2DS4")
nksig_genes <- c("KLRC1", "GNLY", "TRDC", "FGFBP2", "KLRB1", "KLRC2")
cytotoxic_genes <- c("PRF1", "GZMB", "GZMA", "GZMH", "NKG7", "GNLY")
inhib_genes <- c("PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT")
trm_genes <- c("CA10", "ITGA1", "ITGAE", "IL2", "IL10", "CXCR6", "CXCL13", "KCNK5", "RGS1", "CRTAM", "DUSP6", "PDCD1")
cd8_genes <- c("CD8A", "KLRK1", "NKG7", "CD8B", "KLRC4", "KLRD1", "ZNF683")
cd4_genes <- c("CD4", "CD40LG", "KLRB1", "RUNX2", "IL7R")
lym.obj <- AddModuleScore(lym.obj,
                          features = list(nksig_genes),
                          name = "ModuleScore_NKSig",
                          assay = "RNA")
lym.obj <- AddModuleScore(lym.obj,
                          features = list(cytotoxic_genes),
                          name = "ModuleScore_cytotoxic",
                          assay = "RNA")
lym.obj <- AddModuleScore(lym.obj,
                          features = list(cd4_genes),
                          name = "ModuleScore_CD4",
                          assay = "RNA")
FeaturePlot(lym.obj,
            features = c("ModuleScore_CD41")) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu"))) +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        legend.position = "none",
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))

lym.obj <- AddModuleScore(lym.obj,
                          features = list(cd8_genes),
                          name = "ModuleScore_CD8",
                          assay = "RNA")
FeaturePlot(lym.obj,
            features = c("ModuleScore_CD81")) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "RdBu"))) +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        legend.position = "none",
        plot.title = element_blank(),
        plot.subtitle = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))

Idents(object = lym.obj) <- "seurat_clusters"
markers <- FindAllMarkers(lym.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top7 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 7, wt = avg_log2FC)

lym.cluterHeatmap <- AverageExpression(lym.obj, 
                                       assays = "integrated",
                                       features = c(top7$gene),
                                       group.by = "seurat_clusters",
                                       return.seurat = TRUE)

col_fun = colorRamp2(c(-1.5, 0, 1.5), c("navy", "white", "red"))
Heatmap(lym.cluterHeatmap@assays$integrated@scale.data,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE, 
                   border = "black",
                   column_split = c(0:(ncol(lym.cluterHeatmap@assays$integrated@scale.data)-1)),
                   col = col_fun,
                   show_column_names = FALSE, 
                   heatmap_legend_param = list(
                     title = "Relative expression",
                     direction = "horizontal",
                     at = c(-2, 0, 2),
                     border = "black",
                     legend_width = unit(6, "cm"),
                     title_position = "topcenter"
                   ))

# Cluster 1 expressed myeloid cells, general T-cell genes, and CD4 and CD8 T-Cells --> remove
Idents(object = lym.obj) <- "seurat_clusters"
lym.obj <- subset(x = lym.obj, idents = c(0,2))
cd8TCell <- colnames(lym.obj)[lym.obj$seurat_clusters == 0]
cd4TCell <- colnames(lym.obj)[lym.obj$seurat_clusters == 2]
# KLRB1 can stratify pediatric LGG patient survival if split by high and low expression
FeaturePlot(lym.obj, features = "KLRB1")


