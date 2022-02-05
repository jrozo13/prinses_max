# ESTIMATE funciton analysis
# Last updated: 23/12/2021

########## Organize data and cohort ##########
library(tidyverse)
pa_meta_data <- read.csv(file = "/Users/jrozowsky/Documents/PMC/PA/PA_DataSets/PA_cohort_metadata.csv")
pa_samples <- pa_meta_data$Sample.ID
covariate_data <- meta_data %>% column_to_rownames(var = "Sample.ID") %>% 
  select("gender", "Age..years.", "location_updated", "Primary.or.recurrent.")
colnames(covariate_data) <- c("Gender", "Age", "Location", "Stage")

load_data <- readRDS("/Users/jrozowsky/Documents/PMC/Data/20211126_PMCdiag_RNAseq_counts_noHiX.rds")
meta_data <- load_data$metaData
gbm_samples <- rownames(meta_data)[meta_data$Disease_sub_class == "Glioblastoma"]
dmg_samples <- rownames(meta_data)[meta_data$Disease_sub_class == "Diffuse midline glioma"]
epn_samples <- rownames(meta_data)[meta_data$Disease_sub_class == "Ependymoma"]
mb_samples <- rownames(meta_data)[meta_data$Disease_sub_class == "Medulloblastoma"]
lgg_samples <- setdiff(rownames(meta_data)[meta_data$Disease_main_class == "Low grade glioma" & meta_data$Disease_sub_class == "Pilocytic astrocytoma"], pa_samples)
meta_samples <- c(pa_samples, gbm_samples, dmg_samples, epn_samples, mb_samples, lgg_samples)

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

neuro_estimate %>% 
  select(ImmuneScore, Diagnosis) %>%
  ggplot(aes(x = Diagnosis, y = ImmuneScore)) +
  geom_violin(aes(fill = Diagnosis), scale = "width", width = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
  scale_fill_manual(name = "Diagnosis", values = col_annotations$Diagnosis) +
  scale_color_manual(name = "Diagnosis", values = col_annotations$Diagnosis) +
  stat_compare_means(label.y = 1400, label.x = 1) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank())
ggsave(filename = "/Users/jrozowsky/Documents/PMC/PA/PA_Figures/ImmuneScore_Diagnosis.png", plot = last_plot(),
       width = 10, height = 7)

pa_estimate <- neuro_estimate[pa_samples,] %>%
  merge(covariate_data, by = "row.names") %>%
  column_to_rownames(var = "Row.names")

pa_estimate %>% 
  select(ImmuneScore, Location) %>%
  ggplot(aes(x = Location, y = ImmuneScore)) +
  geom_violin(aes(fill = Location), scale = "width", width = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
  scale_fill_manual(name = "Location", values = col_annotations$Location) +
  stat_compare_means(label.y = 1400, label.x = 1) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank())
ggsave(filename = "/Users/jrozowsky/Documents/PMC/PA/PA_Figures/ImmuneScore_PA_Location.png", plot = last_plot(),
       width = 10, height = 7)
