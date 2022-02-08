# Pilocytic astrocytoma clinical data analysis
# Last updated: 04/02/2022
# Note: skip to bottom to run OncoPrint

########## Cohort ########## 
load_data <- readRDS("/Users/jrozowsky/Library/Mobile Documents/com~apple~CloudDocs/Documents/PMC/Data/PMC/20211126_PMCdiag_RNAseq_counts_noHiX.rds")
meta_data <- load_data$metaData
lennart_samples <- rownames(meta_data)[meta_data$Disease_sub_class == "Pilocytic astrocytoma"]

library(readxl)
mutation_cohort <- read_excel(path = "/Users/jrozowsky/Library/Mobile Documents/com~apple~CloudDocs/Documents/PMC/PA/Cohort/ClinicalData_updates.xlsx",
           sheet = "MutationData_cohort") %>% data.frame()
mutation_cohort$...1 <- NULL
mat <- mutation_cohort %>% 
  dplyr::filter(Unique_samples == "Keep") %>%
  select(Sample.ID, Gene.mutated, Mutation.class) %>%
  spread(key = Gene.mutated, value = Mutation.class, drop = TRUE) %>%
  `rownames<-`(.[,1]) %>% 
  select(-Sample.ID) %>%
  as.matrix() %>%
  t()

mat[is.na(mat)] <- ""
unique_samples <- colnames(mat)

annotations <- read_excel(path = "/Users/jrozowsky/Library/Mobile Documents/com~apple~CloudDocs/Documents/PMC/PA/Cohort/ClinicalData_updates.xlsx",
                          sheet = "ClinicalData_cohort") %>% data.frame()
annotations$...1 <- NULL
annotations_updated <- annotations %>% 
  dplyr::filter(Sample.ID %in% unique_samples) %>%
  dplyr::filter(Color..Mariette. != "non piloA LGG: both morphology or molecular does not fit") %>%
  select(Sample.ID, gender, Age..years., location, X..Mutations)

length(intersect(annotations_updated$Sample.ID, lennart_samples))
# 88 similar samples between Mariette's and Lennart's samples

# Find which samples are missing from Lennart data set
lennart_missing <- setdiff(unique_samples, lennart_samples)
# 1 sample in Mariette's data but not in Lennart's data
# These are annotated as "Low grade glioma, other" --> add into meta_data manually
lennart_missing %in% row.names(meta_data) # TRUE --> we have the data, just annotated wrong
lennart_samples <- unique(c(lennart_samples, "PMLBM000AGP", "PMLBM000DUQ"))
setdiff(unique_samples, unique(lennart_samples)) # fixed!
unique_samples <- intersect(unique_samples, lennart_samples) # 90 samples

# Find which samples are missing from Mariette's data set
mariette_missing <- setdiff(lennart_samples, unique_samples)
# 13 samples in Lennart's data but not in Mariette's
track_missing <- data.frame(matrix(nrow = length(mariette_missing), ncol = 1))
colnames(track_missing) <- c("fixed?")
rownames(track_missing) <- mariette_missing

track_missing["PMABM000GDY",] <- "Atypical PA"
track_missing["PMABM000GJI",] <- "DLGNT"
track_missing["PMLBM000ACK",] <- "Matched recurrence"
track_missing["PMLBM000AGZ",] <- "Matched recurrence"
track_missing["PMLBM000BML",] <- "Atypical PA"
track_missing["PMABM000GIL",] <- "DLGNT"
track_missing["PMLBM000BUE",] <- "Matched recurrence"
track_missing["PMLBM000CJY",] <- "Matched recurrence"
track_missing["PMLBM000CWJ",] <- "DLGNT"
track_missing["PMLBM000CWL",] <- "Matched recurrence"
track_missing["PMLBM000CZG",] <- "Matched recurrence"
track_missing["PMLBM000DRT",] <- "Atypical PA"
track_missing["PMLBM000DWM",] <- "DLGNT"

table(track_missing$`fixed?`)
# Do we include the 4 Atypical PA?
# Recurrent matched sample will not be included in the initial cohort (unique cases)
unique_samples <- c(unique_samples) %>% unique()
meta_pa <- meta_data[unique_samples,]

########## Organize dataframe for OncoPrint ########## 
# make oncoprint of all unique samples
# include those that are recurrent, but not ones that are matched
# order based on: 1) gene 2) # mutations 3) location
combined <- merge(annotations_updated, mutation_cohort, by = "Sample.ID") 
combined$Gene.mutated <- factor(combined$Gene.mutated, 
                                levels = names(rev(sort(table(combined$Gene.mutated)))))
combined$Mutation.class <- factor(combined$Mutation.class, 
                                  levels = c("Fusion; SNV", "Fusion", "Duplication; SNV",
                                  "Duplication", "Deletion", "SNV"))
combined$X..Mutations <- factor(combined$X..Mutations, levels = c("2", "1"))
combined <- combined %>% 
  mutate(location_updated = ifelse(location == "PF", "Posterior fossa",
                                   ifelse(location == "ST", "Supratentorial",
                                          ifelse(location %in% c("craniospinal", "intramedullary"), "Spinal", NA)))) %>%
  mutate(Age..years. = as.numeric(Age..years.)) %>%
  mutate(age_updated = ifelse(Age..years. > 10, "Over 10 years", "Under 10 years"))

combined$location_updated <- factor(combined$location_updated, 
                                  levels = c("Posterior fossa", "Supratentorial", "Spinal", NA))
combined <- combined %>% 
  arrange(Gene.mutated, X..Mutations, Mutation.class, location_updated, Age..years., gender) %>%
  distinct(Sample.ID, .keep_all = TRUE)
#write.csv(combined, file = "/Users/jrozowsky/Library/Mobile Documents/com~apple~CloudDocs/Documents/PMC/PA/PA_Data/PA_cohort_metadata.csv")

sample_order <- unique(combined$Sample.ID)
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
color_set <- brewer.pal(4, "Set1")
col_annotations <- list(
  Sex = c("M" = "skyblue", "F" = "pink"),
  Location = c("Spinal" = "#f2f0f7", 
               "Posterior fossa" = "#6a51a3", 
               "Supratentorial" = "#9e9ac8"),
  `Age (years)` = colorRamp2(c(0,20), c("#edf8e9", "#74c476"))
)

mat <- read.csv("PA_DataSets/PA_for_oncoprint.csv", header = T, row.names = 1)
combined <- read.csv("/Users/jrozowsky/Documents/PMC/PA/PA_DataSets/PA_cohort_metadata.csv",
                     header = T, row.names = 1)

age <- combined$Age..years.
age_updated <- combined$age_updated
location <- combined$location_updated
sex <- combined$gender

# save figure as 5 x 10
oncoprint_PA <- oncoPrint(mat,
          alter_fun = alter_fun, 
          col = col,
          pct_side = "right", 
          row_names_side = "left",
          alter_fun_is_vectorized = FALSE,
          column_order = sample_order,
          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                             `Age (years)` = age,
                                             Location = location,
                                             Sex = sex,
                                             col = col_annotations))
