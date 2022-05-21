# Single-cell Pilocytic Astrocytoma analysis: make Taylor datasets
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

########## Load Rozowsky data ##########
library(Seurat)
library(SingleCellExperiment)

taylor_ids <- list.files(path = paste0(wd, "Data/Vladoiu_2019/"))
# loop through all Taylor files and make seurat object for each sample
for (patient in taylor_ids) {
  file_dir <- paste0(wd, "Data/Vladoiu_2019/", patient, "/")
  counts <- Read10X(data.dir = file_dir)
  seuratObj <- CreateSeuratObject(counts = counts, project = patient, assay = "RNA")
  
  location <- "PF"
  seuratObj$updated_location <- location
  
  print(paste0(patient, "_sce"))
  assign(paste0(patient, "_sce"), seuratObj)
  
  rm(seuratObj)
  rm(counts)
}

taylor.combined <- merge(get(paste0(taylor_ids[1], "_sce")),
                         y = mget(paste0(taylor_ids[-1], "_sce")),
                         add.cell.ids = c(taylor_ids),
                         project = "Taylor_PiloAstro")
rm(list = paste0(taylor_ids, "_sce"))

########## Quality Control ########## 
library(scater)
# change this to the object of interest
# seurat.object = combined
seurat.object = taylor.combined

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
# filter out 0 nuclei

# Quality check: number of genes expressed (nFeature_RNA)
summary(sce$nFeature_RNA)
ggplot(mapping = aes(x = sce$nFeature_RNA)) + 
  geom_density(fill = "lightblue") +
  geom_vline(xintercept = 500, color = "red") +
  labs(x = "Number of genes expressed")
sum(sce$nFeature_RNA < 500)
# Filter out 16 nuclei

# Qualtiy check: percent mitochondrial genes expressed (percent.mito)
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
sum(sce$percent.mito > 20)
# Filter 63 nuclei

# Set Thresholds for qualtiy control
genes_exp.check = 500
total_counts.check = 1000
mito_fraction.check = 20

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

assign("Taylor_sce", sce_clean)

save(Taylor_sce, file = paste0(wd, "PA/PA_Data/Taylor.SingleCellExpObject_19.05.2022.RData"))

