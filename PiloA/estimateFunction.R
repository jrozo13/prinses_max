# ESTIMATE funciton analysis
# Last updated: 23/12/2021

########## Organize data and cohort ##########
library(tidyverse)
pa.meta_data <- read_excel(path = "PA/Cohort/ClinicalData_updates.xlsx", sheet = "CohortList")
pa_samples <- pa.meta_data$PMABM

load_data <- readRDS("/Users/jrozowsky/Documents/PMC/Data/PMC/20211126_PMCdiag_RNAseq_counts_noHiX.rds")
meta_data <- load_data$metaData
gbm_samples <- rownames(meta_data)[meta_data$Disease_sub_class == "Glioblastoma"]
dmg_samples <- rownames(meta_data)[meta_data$Disease_sub_class == "Diffuse midline glioma"]
epn_samples <- rownames(meta_data)[meta_data$Disease_sub_class == "Ependymoma"]
mb_samples <- rownames(meta_data)[meta_data$Disease_sub_class == "Medulloblastoma"]
lgg_samples <- setdiff(rownames(meta_data)[meta_data$Disease_main_class == "Low grade glioma" & meta_data$Disease_sub_class == "Pilocytic astrocytoma"], pa_samples)
meta_samples <- c(pa_samples, gbm_samples, dmg_samples, epn_samples, mb_samples, lgg_samples)
meta_samples <- meta_samples[meta_samples %in% rownames(meta_data)]

meta_df <- meta_data %>% 
  rownames_to_column() %>%
  select(rowname, Disease_sub_class, Disease_main_class) %>%
  dplyr::filter(rowname %in% meta_samples) %>%
  mutate(Diagnosis = ifelse(rowname %in% pa_samples, "Pilocytic astrocytoma",
                            ifelse(rowname %in% lgg_samples, "Non-PA LGG",
                                   ifelse(rowname %in% gbm_samples, "Glioblastoma",
                                          ifelse(rowname %in% dmg_samples, "Diffuse midline glioma",
                                                 ifelse(rowname %in% epn_samples, "Ependymoma", "Medulloblastoma")))))) %>%
  column_to_rownames(var = "rowname") %>%
  select(Diagnosis)


count_data <- load_data$rawCounts %>% select(meta_samples)
rm(load_data)

########## ESTIMATE ##########
library(estimate)

write.table(data.frame(count_data), file = "/Library/Frameworks/R.framework/Versions/4.1/Resources/library/estimate/extdata/neuro_input.txt",
            sep = "\t", row.names = TRUE, quote = FALSE)
neuro_input <- system.file("extdata", "neuro_input.txt", package="estimate")
neuro_genes <- "/Library/Frameworks/R.framework/Versions/4.1/Resources/library/estimate/extdata/neuro_genes.gct"
filterCommonGenes(input.f = neuro_input, output.f = neuro_genes, id="GeneSymbol")
neuro_output <- "/Library/Frameworks/R.framework/Versions/4.1/Resources/library/estimate/extdata/neuro_estimate.gct"
estimateScore(neuro_genes, neuro_output, platform = "illumina")

neuro_estimate <- read.delim(file = neuro_output, skip = 2, row.names = 2) %>% t() %>% data.frame()
neuro_estimate <- neuro_estimate[-1,]
neuro_estimate$StromalScore <- as.numeric(neuro_estimate$StromalScore)
neuro_estimate$ImmuneScore <- as.numeric(neuro_estimate$ImmuneScore)
neuro_estimate$ESTIMATEScore <- as.numeric(neuro_estimate$ESTIMATEScore)

neuro_estimate <- merge(neuro_estimate, meta_df, by = "row.names") %>%
  column_to_rownames(var = "Row.names")
neuro_estimate$Diagnosis <- factor(neuro_estimate$Diagnosis, levels = 
                              c("Medulloblastoma", "Ependymoma", "Diffuse midline glioma", "Glioblastoma", "Non-PA LGG", "Pilocytic astrocytoma"))

library(ggpubr)
library(circlize)
library(scales)
col_annotations <- list(
  Sex = c("M" = "skyblue", "F" = "pink"),
  Location = c("Spinal" = "#f2f0f7", 
               "Posterior fossa" = "#6a51a3", 
               "Supratentorial" = "#9e9ac8"),
  Age = colorRamp2(c(0,20), c("#edf8e9", "#74c476"),),
  Diagnosis = c("Medulloblastoma" = hue_pal()(6)[1],
                "Ependymoma" = hue_pal()(6)[2],
                "Diffuse midline glioma" = hue_pal()(6)[3],
                "Glioblastoma" = hue_pal()(6)[4],
                "Non-PA LGG" = hue_pal()(6)[5],
                "Pilocytic astrocytoma" = "#6a51a3")
)
library(viridis)
colors = plasma(n = 6)
#Bulk.ImmuneScoreEnrichmentDiagnosis <- 
neuro_estimate %>% 
  select(StromalScore, Diagnosis) %>%
  ggplot(aes(x = Diagnosis, y = StromalScore)) +
  geom_violin(aes(fill = Diagnosis), scale = "width", width = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.7) +
  scale_fill_manual(values = colors) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 15)) +
  coord_flip()
pdf(file = paste0(fwd, "Bulk.ImmuneScoreEnrichmentDiagnosis.pdf"), width = 6, height = 4)
print(Bulk.ImmuneScoreEnrichmentDiagnosis)
dev.off()

pa_estimate <- neuro_estimate[pa_samples,] %>%
  merge(covariate_data, by = "row.names") %>%
  column_to_rownames(var = "Row.names")

pa_estimate %>% 
  select(ImmuneScore, Location) %>%
  ggplot(aes(x = Location, y = ImmuneScore)) +
  geom_violin(aes(fill = Location), scale = "width", width = 0.7) +
  geom_boxplot(width = 0.2, fill = "white", alpha = 0.7) +
  scale_fill_manual(values = colors) +
  stat_compare_means(label.y = 1400, label.x = 1) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank())
ggsave(filename = "/Users/jrozowsky/Documents/PMC/PA/PA_Figures/ImmuneScore_PA_Location.png", plot = last_plot(),
       width = 10, height = 7)
