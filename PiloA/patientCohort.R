# Pilocytic astrocytoma clinical data analysis
# Last updated: 04/02/2022
# Note: skip to bottom to run OncoPrint

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
  Location = c("spinal" = "#EF553B", 
               "posterior fossa" = "#636EFA", 
               "supratentorial" = "#00CC96"),
  Age = colorRamp2(c(0,20), c("#edf8e9", "#74c476"))
)

wd <- "/Users/jrozowsky/Library/Mobile Documents/com~apple~CloudDocs/Documents/PMC/"
setwd(wd)
fwd <- paste0(wd, "PA/Analysis/Figures/") # figure working directory


########## Cohort ########## 
library(readxl)
cohort <- read_excel(path = "/Users/jrozowsky/Library/Mobile Documents/com~apple~CloudDocs/Documents/PMC/PA/Cohort/ClinicalData_updates.xlsx",
                sheet = "CohortList") %>% data.frame()
mat <- cohort %>%
  filter(UniqueSample == "Yes") %>%
  select(Subject.ID, Molecular.1, Molecular.2) %>%
  separate_rows(Molecular.1, Molecular.2, sep = ",") %>%
# when there are multiple mutated genes, separate by comma in excel sheet
# when there are muttiple mutations in one gene, separate mutation types by colon
# this way, can keep the colon separator in the final mutation matrix
  spread(key = Molecular.1, value = Molecular.2, drop = TRUE) %>% 
  column_to_rownames(var = "Subject.ID") %>%
  as.matrix() %>%
  t()

mat[is.na(mat)] <- ""
unique_samples <- colnames(mat)

annotations <- cohort %>%
  filter(UniqueSample == "Yes") %>%
  select(Subject.ID, Sex, AgeMonths, Status, Location, LocationSpecific, Molecular.1, Molecular.2, NumMut, snRNAseq)

########## Organize dataframe for OncoPrint ########## 
# make oncoprint of all unique samples
# include those that are recurrent, but not ones that are matched
# order based on: 1) gene 2) # mutations 3) location
annotations <- annotations %>%
  mutate(Molecular.2 = factor(annotations$Molecular.2,
                              levels = c("Fusion",
                                         "Deletion",
                                         "Duplication",
                                         "SNV",
                                         "Fusion, SNV",
                                         "Fusion; SNV",
                                         "SNV, SNV"))) %>%
  mutate(Location = factor(annotations$Location,
                           levels = c("PF", "ST", "intramedullary"))) %>%
  mutate(snRNAseq = ifelse(snRNAseq == "No", "No", "Yes"))
  

library(plyr)
annotations$Location <- revalue(annotations$Location,
                                c("intramedullary" = "spinal",
                                  "PF" = "posterior fossa",
                                  "ST" = "supratentorial"))
annotations <- annotations %>% arrange(Molecular.1, Molecular.2, NumMut, Location, Sex)

#write.csv(combined, file = "/Users/jrozowsky/Library/Mobile Documents/com~apple~CloudDocs/Documents/PMC/PA/PA_Data/PA_cohort_metadata.csv")

sample_order <- unique(annotations$Subject.ID)
mat <- mat[,sample_order]
#write.csv(mat, file = "/Users/jrozowsky/Library/Mobile Documents/com~apple~CloudDocs/Documents/PMC/PA/PA_Data/PA_for_oncoprint.csv")

########## Make OncoPrint ##########
library(xlsx)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
color_set <- brewer.pal(4, "Set1")
col = c(SNV = color_set[1], 
        Fusion = color_set[2],
        Duplication = color_set[3],
        Deletion = color_set[4])
alter_fun = list(
  background = alter_graphic("rect", fill = "#CCCCCC"), 
  Duplication = alter_graphic("rect", fill = col["Duplication"]),
  Fusion = alter_graphic("rect", fill = col["Fusion"]),
  Deletion = alter_graphic("rect", fill = col["Deletion"]),
  SNV = alter_graphic("rect", height = 0.33, fill = col["SNV"]))

# mat <- read.csv("/Users/jrozowsky/Library/Mobile Documents/com~apple~CloudDocs/Documents/PMC/PA/PA_Data/PA_for_oncoprint.csv", header = T, row.names = 1)
combined <- read.csv(paste0(wd, "PA/PA_Data/PA_cohort_metadata.csv"),
                     header = T, row.names = 1)

age <- annotations$AgeMonths
location <- annotations$Location
location.2 <- annotations$LocationSpecific
sex <- annotations$Sex

RNAseq <- rep(c("Yes"), 83)
DNAme <- rep(c("Yes"), 83)
snRNAseq <- annotations$snRNAseq
bottom_mat <- data.frame(RNAseq, DNAme, snRNAseq)
bottom_mat[bottom_mat == "Yes"] <- 20
bottom_mat[bottom_mat == "No"] <- NA
bottom_col = c("Yes" = "white", "No" = "white")

# save figure as 5 x 10
oncoprint_PA <- oncoPrint(mat,
          alter_fun = alter_fun, 
          col = col,
          pct_side = "right", 
          pct_gp = gpar(fontsize = 15),
          row_names_side = "left",
          alter_fun_is_vectorized = FALSE,
          column_order = sample_order,
          row_names_gp = gpar(fontsize=15),
          top_annotation = HeatmapAnnotation(
            cbar = anno_oncoprint_barplot(),
            #Age = age,
            Location = location,
            Sex = sex,
            col = col_annotations,
            annotation_name_gp = gpar(fontsize = 15)),
          bottom_annotation = HeatmapAnnotation(
            DNAme = anno_simple(RNAseq,
                                pch = as.numeric(bottom_mat$DNAme),
                                pt_gp = gpar(col = "black"),
                                pt_size = unit(rep(c(2),83),"mm"),
                                col = bottom_col),
            RNAseq = anno_simple(DNAme,
                                 pch = as.numeric(bottom_mat$RNAseq),
                                 pt_gp = gpar(col = "black"),
                                 pt_size = unit(rep(c(2),83),"mm"),
                                 col = bottom_col),
            snRNAseq = anno_simple(snRNAseq,
                                 pch = as.numeric(bottom_mat$snRNAseq), 
                                 pt_gp = gpar(col = "black"), 
                                 pt_size = unit(rep(c(2),83),"mm"),
                                 col = bottom_col),
            annotation_height = unit(10, "mm"),
            annotation_name_gp = gpar(fontsize = 10)))
            
pdf(paste0(fwd, "oncoplot_22.05.22.pdf"), width = 14, height = 5.5)
print(oncoprint_PA)
dev.off()
