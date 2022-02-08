# Pilocytic astrocytoma transcriptomic analysis
# Last updated: 11/01/2022

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
  Location = c("Spinal" = "#EF553B", 
               "Posterior fossa" = "#636EFA", 
               "Supratentorial" = "#00CC96"),
  Age = colorRamp2(c(0,20), c("#edf8e9", "#74c476"))
)

meta_data <- read.csv(file = "/Users/jrozowsky/Library/Mobile Documents/com~apple~CloudDocs/Documents/PMC/PA/PA_Data/PA_cohort_metadata.csv")
unique_samples <- meta_data$Sample.ID
covariate_data <- meta_data %>% column_to_rownames(var = "Sample.ID") %>% 
  select("gender", "Age..years.", "location_updated", "Primary.or.recurrent.") %>%
  drop_na()
colnames(covariate_data) <- c("Sex", "Age", "Location", "Stage")
covariate_data$Location <- as.factor(covariate_data$Location)
covariate_data$Sex <- as.factor(covariate_data$Sex)
unique_samples <- rownames(covariate_data)

load_data <- readRDS("/Users/jrozowsky/Library/Mobile Documents/com~apple~CloudDocs/Documents/PMC/Data/PMC/20211126_PMCdiag_RNAseq_counts_noHiX.rds")
count_data <- load_data$rawCounts %>% select(unique_samples)
cpm_data <- load_data$counts %>% data.frame() %>% select(unique_samples)
rm(load_data)

########## Whole transcriptome analysis ##########
countsLog <- log(cpm_data + 1)
varGenes <- apply(countsLog, 1, var)
meanGenes <- apply(countsLog, 1, mean)
nFeatures = 5000
varFeatures <- names(varGenes)[order(varGenes, decreasing = T)][c(1:nFeatures)]
dataScale <- apply(countsLog, 2, function(x) (x - meanGenes)/varGenes)

### Dimensionality reduction ###
### PCA ###
library(DESeq2)
# library("gg3D")
library(plotly)
dds_pca <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = covariate_data,
                              design = ~ 1)
levels(dds_pca$Location) <- c("Posterior.fossa", "Spinal", "Supratentorial")
design(dds_pca) <- ~Location
dds <- DESeq(dds_pca, test = "LRT", reduced = ~ 1)
vsd <- vst(dds_pca, blind = FALSE)
pca <- prcomp(t(assay(vsd)[varFeatures[1:500], ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
d <- cbind(pca$x, covariate_data)

plot_ly(d, x = d$PC1, y = d$PC3, z = d$PC4, color = d$Location, colors = col_annotations$Location) %>%
  add_markers()

ggplot(data = d, aes(x = PC1, y = PC3, color = Location)) +
  geom_point()
ggplot(data = d, aes(x = PC1, y = PC3, color = Location)) +
  geom_point()
ggplot(data = d, aes(x = PC3, y = PC4, color = Location)) +
  geom_point()

# PC1: separates PF from Spinal and ST
# PC2: separates M from F
# PC3: separates ST from PF and Spinal
# PC4: separates Spinal from ST and PF

### UMAP ###
library(umap)
library(ggforce)
set.seed(14); umap_res <- umap(pca$x, alpha=0.1, gamma=0.5)
umap_data <- cbind(umap_res$layout, covariate_data)
ggplot(data = umap_data, aes(x = `1`, y = `2`, color = Location)) +
  #geom_density2d(alpha = 0.7, na.rm = TRUE, contour_var = "count") +
  geom_point(aes(fill = Location), colour = "grey30", pch = 21, size = 2) + 
  xlab("UMAP1") + ylab("UMAP2") +
  guides(col = guide_legend(ncol = 1)) +
  scale_color_manual(name = "Location", values = col_annotations$Location) +
  scale_fill_manual(name = "Location", values = col_annotations$Location) +
  theme_classic() +
  theme(panel.background = element_rect(colour = "grey30", size=2),
        axis.line = element_blank(),
        plot.title = element_text(size = 20),
        legend.position = "none") +
  ylim(c(-2.7, 2.7)) +
  xlim(c(-2.7, 2.7))
#ggsave(filename = "/Users/jrozowsky/Documents/PMC/PA/PA_Figures/UMAP_PA.png", plot = last_plot(), width = 5, height = 5)

### Unsupervised clustering ###
library(ComplexHeatmap)
# Identify protein-coding genes
annotations_ahb <- read.csv("/Users/jrozowsky/Library/Mobile Documents/com~apple~CloudDocs/Documents/PMC/Methods/FunctionalEnrichment/annotations_ahb.csv")
coding_genes <- annotations_ahb$gene_name[annotations_ahb$gene_biotype == "protein_coding"]
varFeatures_coding <- varFeatures[varFeatures %in% coding_genes]

# calculate inter-patient correlation based on protein-coding gene expression
corr_res <- cor(dataScale[varFeatures_coding[1:420],], method = "pearson")
Heatmap(corr_res,
        show_heatmap_legend = FALSE,
        clustering_distance_rows = "pearson",
        clustering_distance_columns = "pearson",
        column_dend_height = unit(4, "cm"),
        top_annotation = HeatmapAnnotation(Age = covariate_data$Age,
                                           Sex = covariate_data$Sex,
                                           Location = covariate_data$Location,
                                           col = col_annotations, 
                                           show_legend = TRUE),
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_row_dend = FALSE)

### Consensus clustering ###
library(ConsensusClusterPlus)
d <- count_data
mads=apply(d,1,mad)
d=d[rev(order(mads))[1:5000],]
d = sweep(d,1, apply(d,1,median,na.rm=T)) %>% as.matrix()
results = ConsensusClusterPlus(d,
                               maxK=6,
                               reps=50,
                               pItem=0.8,
                               pFeature=1,
                               title="title",
                               clusterAlg="hc",
                               distance="pearson",
                               seed= 5,
                               plot=NULL)

########## Differential expression analysis ##########
library(DESeq2)
library(fgsea)
library(org.Hs.eg.db)
library(clusterProfiler)
annotations_ahb <- read.csv("/Users/jrozowsky/Documents/PMC/PA/PA_Methods/annotations_ahb.csv")

#### Spinal vs Supratentorial ####
res.sp_st <- results(dds, alpha = 0.05, test = "Wald", contrast = c("Location", "Spinal", "Supratentorial"|"Posterior.fossa"))
res_tab.sp_st <- data.frame(res.sp_st[order(res.sp_st$padj)[1:5000],])
res_tab.sp_st %>% 
  mutate(threshold = padj < 0.001 & abs(log2FoldChange) >= 1) %>%
  rownames_to_column() %>%
  mutate(text = ifelse(threshold == TRUE, rowname, NA)) %>%
  column_to_rownames(var = "rowname") %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = threshold, label = text)) +
  geom_point(alpha = 0.3) +
  geom_text_repel(box.padding = 0.1, size = 2) +
  ggtitle("DE Genes for Spinal vs Supratentorial PAs") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

ranks <- res_tab.spinal %>% 
  rownames_to_column() %>%
  dplyr::select(rowname, stat) %>%
  na.omit() %>% 
  distinct() %>% 
  group_by(rowname) %>% 
  summarize(stat=mean(stat)) %>%
  deframe()

pathways.hallmark <- gmtPathways("/Users/jrozowsky/Documents/PMC/PA/PA_DataSets/h.all.v7.4.symbols.gmt")
fgseaRes <- fgsea(pathways = pathways.hallmark, stats = ranks, nperm = 1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA") + 
  theme_minimal()

#### Spinal vs PF ####
res.sp_pf <- results(dds, alpha = 0.05, test = "Wald", contrast = c("Location", "Spinal", "Posterior.fossa"))
res_tab.sp_pf <- data.frame(res.sp_pf[order(res.sp_pf$padj)[1:5000],])
res_tab.sp_pf %>% 
  mutate(threshold = padj < 0.001 & abs(log2FoldChange) >= 1) %>%
  rownames_to_column() %>%
  mutate(text = ifelse(threshold == TRUE, rowname, NA)) %>%
  column_to_rownames(var = "rowname") %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = threshold, label = text)) +
  geom_point(alpha = 0.3) +
  geom_text_repel(box.padding = 0.1, size = 2) +
  ggtitle("DE Genes for Spinal vs Posterior fossa PAs") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

#### PF vs ST ####
res.st_pf <- results(dds, alpha = 0.05, test = "Wald", contrast = c("Location", "Supratentorial", "Posterior.fossa"))
res_tab.st_pf <- data.frame(res.st_pf[order(res.st_pf$padj)[1:5000],])
res_tab.st_pf %>% 
  mutate(threshold = padj < 0.001 & abs(log2FoldChange) >= 1) %>%
  rownames_to_column() %>%
  mutate(text = ifelse(threshold == TRUE, rowname, NA)) %>%
  column_to_rownames(var = "rowname") %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = threshold, label = text)) +
  geom_point(alpha = 0.3) +
  geom_text_repel(box.padding = 0.1, size = 2) +
  ggtitle("DE Genes for Supratentorial vs Posterior fossa PAs") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))


covariate_data.2 <- covariate_data %>%
  mutate(Spinal_Loc = ifelse(Location == "Spinal", "Spinal", "Not.spinal")) %>%
  mutate(PF_Loc = ifelse(Location == "Posterior fossa", "Posterior fossa", "Not.PF")) %>%
  mutate(ST_Loc = ifelse(Location == "Supratentorial", "Supratentorial", "Not.ST"))
covariate_data.2$Spinal_Loc <- as.factor(covariate_data.2$Spinal_Loc)
covariate_data.2$PF_Loc <- as.factor(covariate_data.2$PF_Loc)
covariate_data.2$ST_Loc <- as.factor(covariate_data.2$ST_Loc)

#### Spinal vs others ####
dds_sp <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = covariate_data.2,
                              design = ~ Spinal_Loc)
dds_sp <- DESeq(dds_sp)
vsd_sp <- vst(dds_sp, blind = FALSE)

res.spinal <- results(dds_sp, alpha = 0.05, test = "Wald", contrast = c("Spinal_Loc", "Spinal", "Not.spinal"))
res_tab.spinal <- data.frame(res.spinal[order(res.spinal$padj)[1:5000],])
res_tab.spinal %>% 
  mutate(threshold = padj < 0.001 & abs(log2FoldChange) >= 1) %>%
  rownames_to_column() %>%
  mutate(text = ifelse(threshold == TRUE, rowname, NA)) %>%
  column_to_rownames(var = "rowname") %>%
  ggplot(aes(x = log2FoldChange, y = -log10(padj), colour = threshold, label = text)) +
  geom_point(alpha = 0.3) +
  geom_text_repel(box.padding = 0.1, size = 2) +
  ggtitle("DE Genes for Spinal vs Non-spinal PAs") +
  xlab("log2 fold change") + 
  ylab("-log10 adjusted p-value") +
  theme(legend.position = "none",
        plot.title = element_text(size = rel(1.5), hjust = 0.5),
        axis.title = element_text(size = rel(1.25)))

res_tab.spinal <- left_join(res_tab.spinal %>% rownames_to_column(), annotations_ahb,  by=c("rowname"="gene_name"))
all_genes.spinal <- as.character(res_tab.spinal$gene_id)

sigOE_spinal <- dplyr::filter(res_tab.spinal, padj < 0.05 & log2FoldChange > 0) %>%
  dplyr::select(gene_id) %>% na.omit()
sigOE_spinal <- as.character(sigOE_spinal$gene_id)

sigUE_spinal <- dplyr::filter(res_tab.spinal, padj < 0.05 & log2FoldChange < 0) %>%
  dplyr::select(gene_id) %>% na.omit()
sigUE_spinal <- as.character(sigUE_spinal$gene_id)

enrichGO(gene = sigOE_spinal, 
         universe = all_genes.spinal,
         OrgDb = org.Hs.eg.db, 
         keyType = 'ENSEMBL',
         ont = "BP", 
         pAdjustMethod = "BH", 
         qvalueCutoff = 0.05, 
         readable = TRUE) %>%
  enrichplot::pairwise_termsim() %>%
  dotplot(showCategory = 20)


#### PF vs others ####
dds_pf <- DESeqDataSetFromMatrix(countData = count_data,
                                 colData = covariate_data.2,
                                 design = ~ PF_Loc)
dds_pf <- DESeq(dds_pf)
vsd_pf <- vst(dds_pf, blind = FALSE)

res.pf <- results(dds_pf, alpha = 0.05, test = "Wald", contrast = c("PF_Loc", "Posterior fossa", "Not.PF"))
res_tab.pf <- data.frame(res.pf[order(res.pf$padj)[1:5000],])
res_tab.pf <- left_join(res_tab.pf %>% rownames_to_column(), annotations_ahb,  by=c("rowname"="gene_name"))
all_genes.pf <- as.character(res_tab.pf$gene_id)

sigOE_pf <- dplyr::filter(res_tab.pf, padj < 0.05 & log2FoldChange > 0) %>%
  dplyr::select(gene_id) %>% na.omit()
sigOE_pf <- as.character(sigOE_pf$gene_id)

sigUE_pf <- dplyr::filter(res_tab.pf, padj < 0.05 & log2FoldChange < 0) %>%
  dplyr::select(gene_id) %>% na.omit()
sigUE_pf <- as.character(sigUE_pf$gene_id)

enrichGO(gene = sigOE_pf, 
         universe = all_genes.pf,
         OrgDb = org.Hs.eg.db, 
         keyType = 'ENSEMBL',
         ont = "BP", 
         pAdjustMethod = "BH", 
         qvalueCutoff = 0.5, 
         readable = TRUE) %>%
  enrichplot::pairwise_termsim() %>%
  dotplot(showCategory = 20)

########## Immune Panel ##########
immune_genes <- read_excel("Documents/PMC/PA/PA_DataSets/LBL-10043-08_nCounter_PanCancer_Immune_Profiling_Panel_Gene_List.xlsx", 
                                                                               sheet = "Annotations", skip = 1)
immune_genes <- immune_genes$`Gene Name`[1:770]
immune_genes <- intersect(immune_genes, varFeatures)

dataScale_immune <- apply(countsLog[immune_genes,], 2, function(x) (x - meanGenes[immune_genes])/varGenes[immune_genes])
remove <- c("ICAM4")
dataScale_immune <- dataScale_immune[!rownames(dataScale_immune) %in% remove, ]

dataLog_immune <- countsLog[immune_genes,] %>% as.matrix()

Heatmap(dataScale_immune,
        clustering_distance_columns = "pearson",
        clustering_distance_rows = "pearson",
        top_annotation = HeatmapAnnotation(Location = covariate_data$Location,
                                           col = col_annotations, 
                                           show_legend = TRUE),
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_row_dend = TRUE)

library(pheatmap)
pheatmap(mat = dataLog_immune,
         clustering_distance_cols = "pearson",
         annotation_col = covariate_data %>% select(Location, Sex),
         annotation_colors = col_annotations,
         scale = "row",
         show_rownames = FALSE,
         show_colnames = FALSE
          )

########## Deconvolution analysis ##########
### Make bulk and single-cell expression sets
scPA_eset <- readRDS("/Users/jrozowsky/Documents/PMC/PA/PA_DataSets/PA_singleCellReitman.RDS")
sc_genes <- rownames(scPA_eset@assayData[["exprs"]])
meta_genes <- unique(intersect(sc_genes, rownames(count_data)))

bulkPA_eset <- ExpressionSet(assayData = as.matrix(count_data[meta_genes,]))

### Run MuSiC ###
library(MuSiC)
Estimate.prop <- music_prop(bulk.eset = bulkPA_eset, 
                            sc.eset = scPA_eset, 
                            clusters = "Cluster", 
                            samples = "Tumor", 
                            verbose = FALSE)

deconv_res <- data.frame(Estimate.prop$Est.prop.weighted) %>% arrange(-Tumor)
names(deconv_res)[4] <- "Macrophage"
sample_order <- deconv_res %>% arrange(-Tumor) %>% rownames()

### Plot deconvolution results ###
celltype_colors <- c("#F8766D", "seagreen3", "steelblue2", "#CF78FF")
names(celltype_colors) <- c("Tumor", "Microglia", "T.cell", "Macrophage")
deconv_all <- deconv_res %>%
  mutate(Sample = as.character(rownames(deconv_res))) %>%
  gather(key = Cluster, value = Proportion, -Sample) %>%
  mutate(Cluster = factor(Cluster, levels = c("Macrophage", "T.cell", "Microglia", "Tumor")),
         Sample = factor(Sample, levels = sample_order))

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
        plot.title = element_text(size = 15)) +
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
# ggsave(filename = "/Users/jrozowsky/Documents/PMC/PA/PA_Figures/T.cell_deconv.png", plot = last_plot(),
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
# ggsave(filename = "/Users/jrozowsky/Documents/PMC/PA/PA_Figures/Microglia_deconv.png", plot = last_plot(),
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
# ggsave(filename = "/Users/jrozowsky/Documents/PMC/PA/PA_Figures/Tumor_deconv.png", plot = last_plot(),
#        width = 10, height = 7)

########## CIBERSORTx ##########
write.csv(count_data, file = "/Users/jrozowsky/Documents/PMC/PA/PA_DataSets/Bulk_PA_forDeconv.csv")

PA_cibersortx_LM22 <- read.csv("Documents/PMC/PA/PA_DataSets/PA_cibersortx_LM22.csv")
PA_lm22 <- PA_cibersortx_LM22[,1:23] %>%
  column_to_rownames(var = "Mixture") %>%
  dplyr::select(!Dendritic.cells.resting) %>%
  as.matrix() %>% t()
covariate_lm22 <- covariate_data[intersect(rownames(covariate_data), colnames(PA_lm22)),]

pheatmap(PA_lm22,
        show_heatmap_legend = FALSE,
        column_dend_height = unit(4, "cm"),
        top_annotation = HeatmapAnnotation(Age = covariate_data$Age,
                                           Sex = covariate_data$Sex,
                                           Location = covariate_data$Location,
                                           col = col_annotations, 
                                           show_legend = TRUE),
        show_column_names = FALSE,
        show_row_dend = FALSE)
pheatmap(mat = PA_lm22,
         annotation_col = covariate_data %>% dplyr::select(Location, Sex),
         annotation_colors = col_annotations,
         show_colnames = FALSE)
