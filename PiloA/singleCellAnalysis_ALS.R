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