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
mito_fraction.check = 30

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
# load(paste0(wd, "PA/PA_Data/ALS.SingleCellExpObject_25.04.2022.RData"))

########## Make Seurat Object ##########
counts <- counts(sce_clean)
rownames(counts) <- rownames(sce_clean)
colnames(counts) <- colnames(sce_clean)
allcell.obj <- CreateSeuratObject(counts = counts)
allcell.obj$updated_location <- sce_clean$updated_location

##### Split Seurat Object by Location #####
# location.obj <- SplitObject(allcell.obj, split.by = "updated_location")
# pf.obj <- location.obj$`Posterior fossa`
# st.obj <- location.obj$Supratentorial
# spinal.obj <- location.obj$Spinal
# rm(location.obj, allcell.obj)

# save(pf.obj, st.obj, spinal.obj, file = paste0(wd, "PA/PA_Data/ALS.LocationSeuratObjs_25.04.2022.RData"))

# Normalization and clustering
# load("PA/PA_Data/ALS.LocationSeuratObjs_25.04.2022.RData")
library(Seurat)
# was splitting object by locaiton before --> not anymore
# seuratObj <- merge(x = pf.obj, y = c(st.obj, spinal.obj))
object.list <- SplitObject(allcell.obj, split.by = "orig.ident")
for (i in 1:length(object.list)) {
  print(paste0("index: ", i, "; nCells: ", length(Cells(object.list[[i]]))))
}
# remove samples with < 30 cells after filtering
object.list[4] <- NULL

for (i in 1:length(object.list)) {
  object.list[[i]] <- NormalizeData(object.list[[i]], verbose = FALSE)
  object.list[[i]] <- FindVariableFeatures(object.list[[i]], nfeatures = 2000, verbose = FALSE)
}

# Find anchors
anchors <- FindIntegrationAnchors(object.list = object.list, dims = 1:30)
# Integrate data
integrated <- IntegrateData(anchorset = anchors)

# seuratObj <- NormalizeData(seuratObj, verbose = F)
# seuratObj <- FindVariableFeatures(seuratObj, verbose = F)
integrated <- ScaleData(integrated, verbose = F)
integrated <- RunPCA(integrated, features = VariableFeatures(integrated))
ElbowPlot(object = integrated, ndims = 50)
# use 20 PCs

dims.use <- 20
integrated <- RunUMAP(integrated, dims = 1:dims.use, verbose=F)
# integrated <- RunTSNE(integrated, dims = 1:dims.use, verbose=F)
DimPlot(object = integrated, group.by = "orig.ident", reduction = "umap")
# DimPlot(object = integrated, group.by = "orig.ident", reduction = "tsne")

assign(paste0("integrated_mito", mito_fraction.check), integrated)
# We do not see separation by patient or batches

DimPlot(integrated,
        pt.size = 0.5,
        group.by = "orig.ident",
        reduction = "umap") +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  labs(title = "Pilocytic Astrocytoma",
       subtitle = "Alex's Lemonade Stand, 2022") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        plot.title = element_text(),
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))
filteredCells <- setdiff(Cells(integrated_mito30), Cells(integrated_mito20))
DimPlot(integrated,
        pt.size = 0.5, 
        cells.highlight = filteredCells, 
        sizes.highlight = 0.2,
        reduction = "umap") +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  labs(title = "Pilocytic Astrocytoma",
       subtitle = "Alex's Lemonade Stand, 2022\n Filtered cells at 30%") +
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
integrated <- FindClusters(integrated, resolution = 1, group.singletons = TRUE)
DimPlot(integrated,
        pt.size = 1,
        #split.by = "orig.ident", 
        reduction = "umap") +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  labs(title = "Posterior Fossa Pilocytic Astrocytoma",
       subtitle = "Alex's Lemonade Stand, 2022") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))

markers <- FindAllMarkers(integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)
DoHeatmap(integrated, features = top5$gene) + NoLegend()

library(SingleR)
library(celldex)
ref.data <- celldex::HumanPrimaryCellAtlasData(ensembl = FALSE)
# ensembl == TRUE: will use the ENSEMBL IDs
cluster_predictions <- SingleR(as.SingleCellExperiment(integrated), 
                               ref = ref.data,
                               labels = ref.data$label.main,
                               clusters = integrated$integrated_snn_res.1)
cluster_idents <- data.frame(cluster_predictions@rownames, cluster_predictions@listData$labels)
cell_idents <- merge(x = integrated@meta.data %>% rownames_to_column() %>%  select(rowname, integrated_snn_res.1), 
                     y = cluster_idents,
                     by.x = "integrated_snn_res.1",
                     by.y = "cluster_predictions.rownames")%>%
  column_to_rownames(var = "rowname")
cell_idents <- cell_idents[rownames(integrated@meta.data),]
integrated <- AddMetaData(integrated, metadata = cell_idents$cluster_predictions.listData.labels, col.name = "SingleR_ClustPred")
DimPlot(integrated, reduction = "umap", group.by = "SingleR_ClustPred")
DimPlot(integrated, reduction = "umap", group.by = "integrated_snn_res.1")
DimPlot(integrated, reduction = "umap", group.by = "SingleR_ClustPred", split.by = "updated_location")

########## Annotation of immune cells ##########
load("PA/PA_Data/Pombo.SeuratObj_30.03.2022.RData")
assign("ref.pombo", seuratObj); rm(seuratObj)
library(SingleR)
library(SingleCellExperiment)
library(Seurat)

Idents(integrated) <- "SingleR_ClustPred"
immune.obj <- subset(x = integrated, idents = setdiff(integrated$SingleR_ClustPred %>% unique(), "Astrocyte"))
DefaultAssay(immune.obj) <- "integrated"
immune.obj <- ScaleData(immune.obj, verbose = F)
immune.obj <- RunPCA(immune.obj, features = VariableFeatures(immune.obj))
dims.use <- 30
immune.obj <- RunUMAP(immune.obj, dims = 1:dims.use, min.dist = 0.5, verbose = F)
DimPlot(immune.obj, group.by = "updated_location", reduction = "umap")
DimPlot(immune.obj, group.by = "orig.ident", reduction = "umap")

# Under-cluster cells and annotate using SingleR and Pombo annotations
immune.obj <- FindNeighbors(immune.obj)
immune.obj <- FindClusters(immune.obj, group.singletons = TRUE, resolution = 1)

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
immune.obj <- AddMetaData(immune.obj, metadata = immune.cell_idents$immune.cluster_predictions.listData.labels, col.name = col.name)
DimPlot(immune.obj, reduction = "umap", group.by = col.name, split.by = "updated_location")

als.umap_singleRpomboannot <- DimPlot(immune.obj,
                                      reduction = "umap",
                                      group.by = "SingleR.Pombo",
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
als.umap_singleRoannot <- DimPlot(immune.obj,
                                  reduction = "umap",
                                  group.by = "SingleR.Annotation",
                                  pt.size = 1) +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  labs(title = "Pilocytic Astrocytoma: immune microenvironment",
       subtitle = "SingleR annotation") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        plot.title = element_text(),
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))
als.umap_seuratclusters <- DimPlot(immune.obj,
                                  reduction = "umap",
                                  group.by = "seurat_clusters",
                                  pt.size = 1) +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  labs(title = "Pilocytic Astrocytoma: immune microenvironment") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        plot.title = element_text(),
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))
als.umap_patient <- DimPlot(immune.obj,
                                   reduction = "umap",
                                   group.by = "orig.ident",
                                   pt.size = 1) +
  xlab("UMAP1") + 
  ylab("UMAP2") +
  labs(title = "Pilocytic Astrocytoma: immune microenvironment") +
  theme(panel.background = element_rect(colour = "black", size = 1),
        plot.title = element_text(),
        plot.subtitle = element_text(hjust = 0.5),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title = element_text(size = 10),
        axis.title.x = element_text(hjust = 0),
        axis.title.y = element_text(hjust = 0))

library(gridExtra)
pdf(paste0(fwd, "scUMAP.ALS_immuneMicroenviSingleR.pdf"), width = 20, height = 8)
grid.arrange(als.umap_singleRpomboannot, als.umap_singleRoannot, ncol = 2)
dev.off()

pdf(paste0(fwd, "scUMAP.ALS_immuneMicroenviClustering.pdf"), width = 20, height = 8)
grid.arrange(als.umap_patient, als.umap_seuratclusters, ncol = 2)
dev.off()

# save(immune.obj, file = paste0(wd, "PA/PA_Data/ALS.ImmuneObject_28.04.2022.RData"))


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

markers <- FindAllMarkers(gam.obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 5, wt = avg_log2FC)

gam.cluterHeatmap <- AverageExpression(gam.obj, 
                                   assays = "integrated",
                                   features = top5$gene,
                                   group.by = "seurat_clusters",
                                   return.seurat = TRUE)

col_fun = colorRamp2(c(-2, 0, 2), c("navy", "white", "red"))
Heatmap(gam.cluterHeatmap@assays$integrated@scale.data,
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

# save(gam.obj, gam.cluterHeatmap, file = paste0(wd, "PA/PA_Data/ALS.GAMObject_28.04.2022.RData"))
