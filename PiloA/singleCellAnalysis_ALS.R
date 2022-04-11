# Single-cell Pilocytic Astrocytoma analysis (Alex's Lemonade Stand)
# Last updated: 29/03/2022

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
  geom_vline(xintercept = 750, color = "red") +
  labs(x = "Counts per cell") +
  xlim(0, 30000)
sum(sce$nUMI <= 1000)

summary(sce$nGene)
ggplot(mapping = aes(x = sce$nGene)) + 
  geom_density(fill = "lightblue") +
  geom_vline(xintercept = 500, color = "red") +
  labs(x = "Number of genes expressed")
sum(sce$nGene <= 500)

ggplot(mapping = aes(x = sce$percent.mito)) +
  geom_density(fill = "lightblue") + 
  labs(x = "Mitchondrial fraction") +
  geom_vline(xintercept = 20, color = "red") +
  theme_classic()

qc_df <- data.frame(barcode = Cells(sce),
                    genes_exp = sce$nGene,
                    total_counts = sce$nUMI,
                    mito_fraction = sce$percent.mito)
ggplot(qc_df, aes (x = total_counts,
                   y = genes_exp, 
                   color = mito_fraction)) +
  geom_point(alpha = 0.5) +
  scale_color_viridis_c() + 
  geom_vline(xintercept = 1000, color = "red") +
  geom_hline(yintercept = 500, color = "red") +
  labs(x = "Total Count",
       y = "Number of Genes Expressed",
       color = "Mitochondrial\nFraction") + 
  theme_bw()

filtered_samples <- qc_df %>%
  dplyr::filter(total_counts >= 1000,
                genes_exp >= 500,
                mito_fraction <= 20)

# Remove outliers (101 cells which had lowly expressed genes)
sce_clean <- sce[,rownames(filtered_samples)]
sce_clean@colData <- sce_clean@colData[Cells(sce_clean),]

save(sce_clean, file = paste0(wd, "PA/PA_Data/ALS.SingleCellExpObject_11.04.2022.RData"))
# load(paste0(wd, "PA/PA_Data/ALS.SingleCellExpObject_31.03.2022.RData"))

########## Make Seurat Object ##########
counts <- counts(sce_clean)
rownames(counts) <- rownames(sce_clean)
colnames(counts) <- colnames(sce_clean)
seuratObj <- CreateSeuratObject(counts = counts)

### Split Seurat Object by Location ###
location.obj <- SplitObject(combined, split.by = "updated_location")
pf.obj <- location.obj$`Posterior fossa`
st.obj <- location.obj$Supratentorial
spinal.obj <- location.obj$Spinal
rm(location.obj)

# Normalization and clustering
seuratObj <- NormalizeData(seuratObj, verbose = F)
seuratObj <- FindVariableFeatures(seuratObj, verbose = F)
seuratObj <- ScaleData(seuratObj, verbose = F)
seuratObj <- RunPCA(seuratObj, features = VariableFeatures(seuratObj))
ElbowPlot(object = seuratObj, ndims = 50)
DimHeatmap(seuratObj, dims = 25:30, cells = 5000, balanced = TRUE)
# use 30 PCs

dims.use <- 30
seuratObj <- RunUMAP(seuratObj, dims = 1:dims.use, verbose=F)
DimPlot(object = seuratObj, group.by = "orig.ident")
# see clear separation by patient --> harmony correct data

# Harmony correct for sample because they clustered by sample, not cell type
library(harmony)
theta.use<-1
seuratObj <- RunHarmony(seuratObj, 
                        group.by.vars = "orig.ident",
                        theta = theta.use,  
                        plot_convergence = TRUE)

#### Run UMAP on the harmony-corrected PCA
seuratObj <- RunUMAP(seuratObj,
                     reduction = "harmony", 
                     dims = 1:dims.use, 
                     reduction.name = "umapHarmony",
                     reduction.key = "umapHarmony")
seuratObj <- RunTSNE(seuratObj,
                     reduction = "harmony",
                     dims = 1:dims.use,
                     reduction.name = "tsneHarmony",
                     reduction.key = "tsneHarmony")

DimPlot(seuratObj,
        pt.size = 1,
        group.by = "updated_location",
        reduction = "umapHarmony") +
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

########## Clustering ##########
seuratObj <- FindNeighbors(seuratObj, 
                           dims = 1:dims.use,
                           annoy.metric = "euclidean")
seuratObj <- FindClusters(seuratObj, resolution = 0.1, group.singletons = TRUE)
DimPlot(seuratObj,
        pt.size = 1,
        reduction = "umapHarmony") +
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

