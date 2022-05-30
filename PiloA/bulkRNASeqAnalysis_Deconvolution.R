# Bulk Pilocytic Astrocytoma analysis
# Last updated: 22/05/2022

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

########## Read Single Cell expression set ##########
load("PA/Analysis/Deconvolution/ALS_SingleCell.PA_eset.RData")

########## Make Bulk expression set ##########
# Try with new data
load("PA/PA_Data/CohortBulkRNAseq_04.04.2022.RData")
shared_genes <- unique(intersect(SingleCell.PA_eset@assayData$exprs %>% rownames(), rownames(pa_counts)))
Bulk.PA_eset <- ExpressionSet(assayData = as.matrix(pa_counts[shared_genes,]))
pa.samples <- colnames(pa_counts)

# Try with old data
# load_data <- readRDS("Data/PMC/20211126_PMCdiag_RNAseq_counts_noHiX.rds")
# count_data <- load_data$rawCounts
# pa.samples %in% colnames(count_data)
# Bulk.PA_eset <- ExpressionSet(assayData = as.matrix(count_data[shared_genes,]))

########## Deconvolution analysis ##########
library(MuSiC)
PA.basis = music_basis(ALS_SingleCell.PA_eset, 
                             clusters = "Cluster", 
                             samples = "Sample")
par(mfrow = c(1, 2))
d1 <- dist(t(log(PA.basis$Disgn.mtx + 1e-6)), method = "euclidean")
hc1 <- hclust(d1, method = "complete" )
plot(hc1, cex = 0.6, hang = -1, main = 'Cluster log(Design Matrix)')
d2 <- dist(t(log(PA.basis$M.theta + 1e-8)), method = "euclidean")
hc2 <- hclust(d2, method = "complete")
clusters.type = list(C1 = 'Glioma', 
                     C2 = 'B cell', 
                     C3 = c('prol. TAM', 'TAM'),
                     C4 = c('T cell', 'DC', 'Monocyte'))
cl.type = as.character(ALS_SingleCell.PA_eset$Cluster)
for(cl in 1:length(clusters.type)){
  cl.type[cl.type %in% clusters.type[[cl]]] = names(clusters.type)[cl]
}
pData(ALS_SingleCell.PA_eset)$Cluster = factor(cl.type, levels = c(names(clusters.type), 'CD-Trans', 'Novel1', 'Novel2'))

# 13 selected cell types
s.mouse = unlist(clusters.type)
s.mouse

Estimate.prop.cluster <- music_prop.cluster(bulk.eset = Bulk.PA_eset,
                                            sc.eset = ALS_SingleCell.PA_eset, 
                                            clusters = 'Cluster', 
                                            groups = 'clusterType', 
                                            samples = 'Sample', 
                                            clusters.type = clusters.type)
### Need to include group.markers --> each cluster has a vector of cell-type markers

pData(ALS_SingleCell.PA_eset)$Cluster = factor(cl.type, levels = c(names(clusters.type), 'CD-Trans', 'Novel1', 'Novel2'))

Estimate.prop <- music_prop(bulk.eset = Bulk.PA_eset, 
                            sc.eset = ALS_SingleCell.PA_eset, 
                            clusters = "Cluster", 
                            samples = "Sample", 
                            verbose = TRUE)

deconv_res <- data.frame(Estimate.prop$Est.prop.weighted) %>% arrange(-Glioma)
#names(deconv_res)[4] <- "Macrophage"
sample_order <- deconv_res %>% arrange(-Glioma) %>% rownames()

### Plot deconvolution results ###
#celltype_colors <- c("#F8766D", "seagreen3", "steelblue2", "#CF78FF")
#names(celltype_colors) <- c("Tumor", "Microglia", "T.cell", "Macrophage")
deconv_all <- deconv_res %>%
  mutate(Sample = as.character(rownames(deconv_res))) %>%
  gather(key = Cluster, value = Proportion, -Sample)

ggplot(data = deconv_all,
       aes(x = Sample, y = Proportion, fill = Cluster)) +
  ggtitle(label = "Deconvolution of bulk piloctic astrocytomas") +
  geom_bar(stat = "identity") +
  geom_col(width = 1) +
  theme_classic() +
  theme(legend.title = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(size = 15)) 
+
  scale_fill_manual(name = "Group", values = celltype_colors)

### Deconvolution results by clinical data ###
library(ggpubr)
deconv_byLocation  <- merge(deconv_all, 
                            data.frame(covariate_data %>% select(Location) %>% rownames_to_column()),
                            by.x = "Sample", by.y = "rowname")

cell_comparisons <- list(c("Spinal", "Posterior fossa"), c("Posterior fossa", "Supratentorial"), c("Spinal", "Supratentorial"))

compare_means(Proportion ~ Location,  data = deconv_byLocation[deconv_byLocation$Cluster == "T.cell",])
ggplot(data = deconv_byLocation[deconv_byLocation$Cluster == "T.cell",],
       aes(x = Location, y = Proportion)) +
  geom_violin(aes(fill = Location), scale = "width", width = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
  scale_fill_manual(name = "Location", values = col_annotations$Location) +
  stat_compare_means(comparisons = cell_comparisons) +
  stat_compare_means(label.y = 0.18, label.x = 1.2) +
  theme_classic() +
  theme(axis.line = element_line(size = 1, linetype = "solid"),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14),
        panel.grid.minor = element_line(colour = "grey90")) +
  scale_y_continuous(name="Proportion of T cells", limits=c(0, 0.19), expand = c(0,0))
# ggsave(filename = paste0(fwd, "T.cell_deconv.png"), plot = last_plot(),
#        width = 10, height = 7)

compare_means(Proportion ~ Location,  data = deconv_byLocation[deconv_byLocation$Cluster == "Microglia",])
ggplot(data = deconv_byLocation[deconv_byLocation$Cluster == "Microglia",],
       aes(x = Location, y = Proportion)) +
  geom_violin(aes(fill = Location), scale = "width", width = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
  scale_fill_manual(name = "Location", values = col_annotations$Location) +
  stat_compare_means(comparisons = cell_comparisons) +
  stat_compare_means(label.y = 0.52, label.x = 1.2) +
  theme_classic() +
  theme(axis.line = element_line(size = 1, linetype = "solid"),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14),
        panel.grid.minor = element_line(colour = "grey90")) +
  scale_y_continuous(name="Proportion of Microglia", limits=c(0, 0.55), expand = c(0,0))
# ggsave(filename = paste0(fwd, "Microglia_deconv.png"), plot = last_plot(),
#        width = 10, height = 7)

compare_means(Proportion ~ Location,  data = deconv_byLocation[deconv_byLocation$Cluster == "Tumor",])
ggplot(data = deconv_byLocation[deconv_byLocation$Cluster == "Tumor",],
       aes(x = Location, y = Proportion)) +
  geom_violin(aes(fill = Location), scale = "width", width = 0.7) +
  geom_boxplot(width = 0.1, fill = "white", alpha = 0.7) +
  scale_fill_manual(name = "Location", values = col_annotations$Location) +
  stat_compare_means(comparisons = cell_comparisons, label.y = c(0.47, 0.42, 0.38)) +
  stat_compare_means(label.y = 0.35, label.x = 1.2) +
  theme_classic() +
  theme(axis.line = element_line(size = 1, linetype = "solid"),
        axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 14),
        panel.grid.minor = element_line(colour = "grey90")) +
  scale_y_continuous(name="Proportion of Tumor cells", limits=c(0.3, 1), expand = c(0,0))
# ggsave(filename = paste0(fwd, "Tumor_deconv.png"), plot = last_plot(),
#        width = 10, height = 7)