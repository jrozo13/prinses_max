# Single-cell Pilocytic Astrocytoma analysis: make Rozowsky datasets
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

########## Load ALS data ##########
library(Seurat)
library(SingleCellExperiment)

# compile all patient IDs
meta_data <- read.delim(file = "Data/ALS_2022/SCPCP000002/single_cell_metadata.tsv") %>%
  filter(diagnosis == "Pilocytic astrocytoma")
als_ids <- meta_data %>%
  pull(scpca_sample_id)
library_ids <- meta_data %>%
  pull(scpca_library_id)
patient_dict <- data.frame(als_ids, library_ids)

st_loc <- c("Basal ganglia", "Left mesial temporal", "Supracellar", "Thalamic", "Thalamus")
pf_loc <- c("Cerebellar", "Posterior fossa")
meta_data <- meta_data %>%
  mutate(updated_location = ifelse(tissue_location %in% st_loc, "ST",
                                   ifelse(tissue_location %in% pf_loc, "PF", "Spinal")))


for (patient in als_ids) {
  library <-  patient_dict$library_ids[patient_dict$als_ids==patient]
  object <- readRDS(paste0(wd, "Data/ALS_2022/SCPCP000002/", patient, "/", library, "_filtered.rds"))
  
  # convert from ensmbl ID to gene symbol
  rownames(object) <- object@rowRanges@elementMetadata$gene_symbol
  
  # remove rows where ensmbl ID didnt map to gene symbol
  object <- object[!is.na(rownames(object)),]
  
  seuratObj <- CreateSeuratObject(counts = assay(object), project = patient, assay = "RNA")
  
  location <- meta_data$updated_location[meta_data$scpca_sample_id == patient]
  seuratObj$updated_location <- location
  
  print(paste0(patient, "_sce"))
  assign(paste0(patient, "_sce"), seuratObj)
  
  rm(seuratObj)
  rm(object)
}

als.combined <- merge(get(paste0(als_ids[1], "_sce")),
                       y = mget(paste0(als_ids[-1], "_sce")),
                       add.cell.ids = c(als_ids),
                       project = "ALS_PiloAstro")
rm(list = paste0(als_ids, "_sce"))

########## Quality Control ########## 
library(scater)
# change this to the object of interest
# seurat.object = combined
seurat.object = als.combined

sce <- as.SingleCellExperiment(x = seurat.object)

# some quality results are in the seurat object metadata, but want to calculate on our own
mito.genes <- grep("^MT-", rownames(sce), value = TRUE) # 37 mitochondrial genes for combined; 13 for Taylor/Rozowsky
sce$percent.mito <- (Matrix::colSums(counts(sce)[mito.genes, ])*100)/Matrix::colSums(counts(sce)) # in percentage, not fraction
sce$staticNr <- 1
dim(colData(sce))
colnames(colData(sce))

# Quality check: number of UMIs (nCount_RNA)
summary(sce$nCount_RNA)
ggplot(mapping = aes(x = sce$nCount_RNA)) + 
  geom_density(fill = "lightblue") +
  geom_vline(xintercept = 1000, color = "red") +
  labs(x = "Counts per cell") +
  xlim(0, 30000)
sum(sce$nCount_RNA < 1000)
# filter out 9988 cells

# Quality check: number of genes expressed (nFeature_RNA)
summary(sce$nFeature_RNA)
ggplot(mapping = aes(x = sce$nFeature_RNA)) + 
  geom_density(fill = "lightblue") +
  geom_vline(xintercept = 500, color = "red") +
  labs(x = "Number of genes expressed")
sum(sce$nFeature_RNA < 500)
# Filter out 13566 cells

# Qualtiy check: percent mitochondrial genes expressed (percent.mito)
ggplot(mapping = aes(x = sce$percent.mito)) +
  geom_density(fill = "lightblue") + 
  labs(x = "Mitchondrial fraction") +
  geom_vline(xintercept = 20, color = "red")
sum(sce$percent.mito > 20)
ggplot(mapping = aes(x = sce$percent.mito)) +
  geom_density(fill = "lightblue") + 
  labs(x = "Mitchondrial fraction") +
  geom_vline(xintercept = 25, color = "red")
sum(sce$percent.mito > 25)
ggplot(mapping = aes(x = sce$percent.mito)) +
  geom_density(fill = "lightblue") + 
  labs(x = "Mitchondrial fraction") +
  geom_vline(xintercept = 30, color = "red")
sum(sce$percent.mito > 30)
# Set threshold to 25%; Filter 15762 cells

# Set Thresholds for qualtiy control
genes_exp.check = 500
total_counts.check = 1000
mito_fraction.check = 25

qc_df <- data.frame(barcode = Cells(sce),
                    genes_exp = sce$nFeature_RNA,
                    total_counts = sce$nCount_RNA,
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

# Remove low quality cells based on parameters
sce_clean <- sce[,rownames(filtered_samples)]
sce_clean@colData <- sce_clean@colData[Cells(sce_clean),]

assign("ALS_sce", sce_clean)

save(ALS_sce, file = paste0(wd, "PA/PA_Data/ALS.SingleCellExpObject_19.05.2022.RData"))

