# Pilocytic astrocytoma transcriptomic analysis
# Last updated: 11/03/2022

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

wd <- "/Users/jrozowsky/Library/Mobile Documents/com~apple~CloudDocs/Documents/PMC/"
setwd(wd)
fwd <- paste0(wd, "PA/Analysis/Figures/") # figure working directory

annotations_ahb <- read.csv("PA/Analysis/GSEA/annotations_ahb.csv")
annotations_ahb %>% filter(gene_biotype == "protein_coding") %>% pull(gene_name) -> protein_genes

########## Prepare Data ##########
meta_data <- read.csv(file = "PA/PA_Data/PA_cohort_metadata.csv")
unique_samples <- meta_data$Sample.ID
covariate_data <- meta_data %>% column_to_rownames(var = "Sample.ID") %>% 
  select("gender", "Age..years.", "location_updated", "Primary.or.recurrent.") %>%
  drop_na()
colnames(covariate_data) <- c("Sex", "Age", "Location", "Stage")
covariate_data$Location <- as.factor(covariate_data$Location)
covariate_data$Sex <- as.factor(covariate_data$Sex)
unique_samples <- rownames(covariate_data)

load_data <- readRDS("Data/PMC/20211126_PMCdiag_RNAseq_counts_noHiX.rds")
count_data <- load_data$rawCounts %>% select(unique_samples)
cpm_data <- load_data$counts %>% data.frame() %>% select(unique_samples)
rm(load_data)

xy_genes <- read.delim2("PA/PA_Data/xy_genes.txt", header = TRUE, sep = "\t", dec = ".")

## remove x and y linked genes
genes_to_use <- setdiff(rownames(count_data), xy_genes$Approved.symbol)

countsLog <- log(cpm_data[genes_to_use,] + 1)
varGenes <- apply(countsLog, 1, var)
meanGenes <- apply(countsLog, 1, mean)
nFeatures = 5000
varFeatures <- names(varGenes)[order(varGenes, decreasing = T)][c(1:nFeatures)]
dataScale <- apply(countsLog, 2, function(x) (x - meanGenes)/varGenes)

# save(count_data, dataScale, varFeatures, covariate_data, 
#      file = paste0(wd, "PA/PA_Data/PiloA_Data_09.02.2022.RData"))

########## Whole transcriptome analysis ##########
load(file = "PA/PA_Data/PiloA_Data_09.02.2022.RData")

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
dds_pca <- DESeq(dds_pca, test = "LRT", reduced = ~ 1)
vsd <- vst(dds_pca, blind = FALSE)

# Select the most variably expressed genes -> these are non-x/y genes
pca <- prcomp(t(assay(vsd)[varFeatures[1:420], ]))
percentVar <- pca$sdev^2/sum(pca$sdev^2)
d <- cbind(pca$x, covariate_data)

pca_3d <- plot_ly(d, x = d$PC1, y = d$PC2, z = d$PC3, color = d$Location, colors = col_annotations$Location) %>%
  add_markers()

jpeg("/Users/jrozowsky/Desktop/")
print(pca_3d)
dev.off()

ggplot(data = d, aes(x = PC1, y = PC2, color = Location)) +
  geom_point()
ggplot(data = d, aes(x = PC1, y = PC3, color = Location)) +
  geom_point()
ggplot(data = d, aes(x = PC3, y = PC4, color = Location)) +
  geom_point()

### UMAP ###
library(umap)
library(ggforce)
set.seed(613); umap_res <- umap(pca$x, alpha=0.5, gamma=1)
umap_data <- cbind(umap_res$layout, covariate_data)
umap_plot <- ggplot(data = umap_data, aes(x = `1`, y = `2`)) +
  geom_point(aes(fill = Location), colour = "grey30", pch = 21, size = 2) + 
  xlab("UMAP1") + ylab("UMAP2") +
  guides(col = guide_legend(ncol = 1)) +
  scale_color_manual(name = "Location", values = col_annotations$Location) +
  scale_fill_manual(name = "Location", values = col_annotations$Location) +
  theme_classic() +
  theme(panel.background = element_rect(colour = "grey30", size=1),
        axis.line = element_blank(),
        plot.title = element_text(size = 20)) +
  ylim(c(-2.7, 2.7)) +
  xlim(c(-2.7, 2.7))

pdf(paste0(fwd, "UMAP.pdf"), width = 8, height = 6)
print(umap_plot)
dev.off()

### Unsupervised clustering ###
library(ComplexHeatmap)
# Identify protein-coding genes
varFeatures_coding <- varFeatures[varFeatures %in% protein_genes]

# calculate inter-patient correlation based on protein-coding gene expression
corr_res <- cor(dataScale[varFeatures_coding[1:510],], method = "pearson")
heatmap <- Heatmap(corr_res,
        show_heatmap_legend = TRUE,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        column_dend_height = unit(4, "cm"),
        row_names_gp = gpar(fontsize = 8),
        top_annotation = HeatmapAnnotation(Age = covariate_data$Age,
                                           Sex = covariate_data$Sex,
                                           Location = covariate_data$Location,
                                           col = col_annotations, 
                                           show_legend = TRUE),
        show_row_names = TRUE,
        show_column_names = FALSE,
        show_row_dend = FALSE)
pdf(paste0(fwd, "heatmap.pdf"))
print(heatmap)
dev.off()

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
## Spinal vs Supratentorial ##
res.sp_st <- results(dds_pca, alpha = 0.05, test = "Wald", contrast = c("Location", "Spinal", "Supratentorial"))
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

### Do not run
ranks <- res_tab.spinal %>% 
  rownames_to_column() %>%
  dplyr::select(rowname, stat) %>%
  na.omit() %>% 
  distinct() %>% 
  group_by(rowname) %>% 
  summarize(stat=mean(stat)) %>%
  deframe()

pathways.hallmark <- gmtPathways("PA/PA_DataSets/h.all.v7.4.symbols.gmt")
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

## Spinal vs PF ##
res.sp_pf <- results(dds_pca, alpha = 0.05, test = "Wald", contrast = c("Location", "Spinal", "Posterior.fossa"))
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

## PF vs Supratentorial ##
res.st_pf <- results(dds_pca, alpha = 0.05, test = "Wald", contrast = c("Location", "Supratentorial", "Posterior.fossa"))
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

########## Functional enrichment analysis ##########
library(fgsea)
library(org.Hs.eg.db)
library(fgsea)
library(clusterProfiler)
library(DESeq2)
load(file = paste0(wd, "PA/PA_Data/PiloA_GSEA_13.02.2022.RData"))
annotations_ahb %>% filter(gene_biotype == "protein_coding") %>% pull(gene_name) -> protein_genes

covariate_data <- covariate_data %>%
  mutate(Spinal_Loc = ifelse(Location == "Spinal", "Spinal", "Not.spinal")) %>%
  mutate(PF_Loc = ifelse(Location == "Posterior fossa", "Posterior fossa", "Not.PF")) %>%
  mutate(ST_Loc = ifelse(Location == "Supratentorial", "Supratentorial", "Not.ST"))
covariate_data$Spinal_Loc <- as.factor(covariate_data$Spinal_Loc)
covariate_data$PF_Loc <- as.factor(covariate_data$PF_Loc)
covariate_data$ST_Loc <- as.factor(covariate_data$ST_Loc)

## Spinal vs others ##
dds_sp <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = covariate_data,
                              design = ~ Spinal_Loc)
dds_sp <- DESeq(dds_sp)
vsd_sp <- vst(dds_sp, blind = FALSE)

res.spinal <- results(dds_sp, alpha = 0.05, test = "Wald", contrast = c("Spinal_Loc", "Spinal", "Not.spinal"))
res_tab.spinal <- data.frame(res.spinal[order(res.spinal$padj)[1:5000],])
res_tab.spinal <- left_join(res_tab.spinal %>% rownames_to_column(), annotations_ahb,  by=c("rowname"="gene_name"))
all_genes.spinal <- as.character(res_tab.spinal$gene_id)

sig_num.spinal <- res_tab.spinal %>%
  filter(padj < 0.001) %>%
  pull(rowname) %>% length()
# 539 genes are differentially expressed in spinal tumors vs other PAs
top15_genes.spinal <- res_tab.spinal %>%
  filter(log2FoldChange > 0) %>%
  filter(rowname %in% protein_genes) %>%
  arrange(padj) %>%
  top_n(n = -15, padj) %>% pull(rowname)

sigOE_spinal <- dplyr::filter(res_tab.spinal, padj < 0.05 & log2FoldChange > 0) %>%
  dplyr::select(gene_id) %>% na.omit()
sigOE_spinal <- as.character(sigOE_spinal$gene_id)
sigUE_spinal <- dplyr::filter(res_tab.spinal, padj < 0.05 & log2FoldChange < 0) %>%
  dplyr::select(gene_id) %>% na.omit()
sigUE_spinal <- as.character(sigUE_spinal$gene_id)

res_tab.spinal %>% 
  mutate(threshold = padj < 0.001) %>%
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

spUE_GO <- enrichGO(gene = sigUE_spinal, 
                    universe = all_genes.spinal,
                    OrgDb = org.Hs.eg.db, 
                    keyType = 'ENSEMBL',
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE) %>%
  enrichplot::pairwise_termsim() %>%
  dotplot(x = "GeneRatio",
          showCategory = 20,
          font.size = 8,
          title = "Pathways under-enriched in spinal tumors")
spOE_GO <- enrichGO(gene = sigOE_spinal, 
                    universe = all_genes.spinal,
                    OrgDb = org.Hs.eg.db, 
                    keyType = 'ENSEMBL',
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE) %>%
  enrichplot::pairwise_termsim() %>%
  dotplot(x = "GeneRatio",
          showCategory = 20,
          font.size = 8,
          title = "Pathways enriched in spinal tumors")
pdf(paste0(fwd, "spinal_GO.pdf"))
par(mfrow = c(1, 2))
print(spUE_GO)
print(spOE_GO)
dev.off()

spinal_gsea_file <- res_tab.spinal %>%
  filter(padj < 0.01) %>%
  select(rowname, log2FoldChange) %>%
  arrange(-log2FoldChange)
write.table(x = spinal_gsea_file, 
            file = paste0(wd, "PA/PA_Data/spinal_geneList.txt"),
            sep = "\t",
            row.names = FALSE)
spinal_gProf_list <- res_tab.spinal %>%
  filter(padj < 0.05 & log2FoldChange >= 0) %>%
  pull(rowname)
spinal_gProf_list[spinal_gProf_list %in% protein_genes] -> spinal_gProf_list
write(x = spinal_gProf_list,
      file = paste0(wd, "PA/PA_Data/spinal_OE_geneList.txt"),
      sep = "\t")

## PF vs others ##
dds_pf <- DESeqDataSetFromMatrix(countData = count_data,
                                 colData = covariate_data,
                                 design = ~ PF_Loc)
dds_pf <- DESeq(dds_pf)
vsd_pf <- vst(dds_pf, blind = FALSE)

res.pf <- results(dds_pf, alpha = 0.05, test = "Wald", contrast = c("PF_Loc", "Posterior fossa", "Not.PF"))
res_tab.pf <- data.frame(res.pf[order(res.pf$padj)[1:5000],])
res_tab.pf <- left_join(res_tab.pf %>% rownames_to_column(), annotations_ahb,  by=c("rowname"="gene_name"))
all_genes.pf <- as.character(res_tab.pf$gene_id)

sig_num.pf <- res_tab.pf %>%
  filter(padj < 0.001) %>%
  pull(rowname) %>% length()
# 1736 genes are differentially expressed in PF tumors vs other PAs
top15_genes.pf <- res_tab.pf %>%
  filter(log2FoldChange > 0) %>%
  filter(rowname %in% protein_genes) %>%
  arrange(padj) %>%
  top_n(n = -15, padj) %>% pull(rowname)

sigOE_pf <- dplyr::filter(res_tab.pf, padj < 0.05 & log2FoldChange > 0) %>%
  dplyr::select(gene_id) %>% na.omit()
sigOE_pf <- as.character(sigOE_pf$gene_id)

sigUE_pf <- dplyr::filter(res_tab.pf, padj < 0.05 & log2FoldChange < 0) %>%
  dplyr::select(gene_id) %>% na.omit()
sigUE_pf <- as.character(sigUE_pf$gene_id)

pfUE_GO <- enrichGO(gene = sigUE_pf, 
                    universe = all_genes.pf,
                    OrgDb = org.Hs.eg.db, 
                    keyType = 'ENSEMBL',
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE) %>%
  enrichplot::pairwise_termsim() %>%
  dotplot(x = "GeneRatio",
          showCategory = 20,
          font.size = 8,
          title = "Pathways under-enriched in posterior fossa tumors")
pfOE_GO <- enrichGO(gene = sigOE_pf, 
                    universe = all_genes.pf,
                    OrgDb = org.Hs.eg.db, 
                    keyType = 'ENSEMBL',
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE) %>%
  enrichplot::pairwise_termsim() %>%
  dotplot(x = "GeneRatio",
          showCategory = 20,
          font.size = 8,
          title = "Pathways enriched in posterior fossa tumors")

pdf(paste0(fwd, "pf_GO.pdf"))
par(mfrow = c(1, 2))
print(pfUE_GO)
print(pfOE_GO)
dev.off()

pf_gsea_file <- res_tab.pf %>%
  filter(padj < 0.01) %>%
  select(rowname, log2FoldChange) %>%
  arrange(-log2FoldChange)
write.table(x = pf_gsea_file, 
            file = "PA/PA_Data/pf_geneList.txt",
            sep = "\t",
            row.names = FALSE)

pf_gProf_list <- res_tab.pf %>%
  filter(padj < 0.05 & log2FoldChange >= 0) %>%
  pull(rowname)
pf_gProf_list[pf_gProf_list %in% protein_genes] -> pf_gProf_list
write(x = pf_gProf_list,
      file = paste0(wd, "PA/PA_Data/pf_OE_geneList.txt"),
      sep = "\t")

## ST vs others ##
dds_st <- DESeqDataSetFromMatrix(countData = count_data,
                                 colData = covariate_data,
                                 design = ~ ST_Loc)
dds_st <- DESeq(dds_st)
vsd_st <- vst(dds_st, blind = FALSE)

res.st <- results(dds_st, alpha = 0.05, test = "Wald", contrast = c("ST_Loc", "Supratentorial", "Not.ST"))
res_tab.st <- data.frame(res.st[order(res.st$padj)[1:5000],])
res_tab.st <- left_join(res_tab.st %>% rownames_to_column(), annotations_ahb,  by=c("rowname"="gene_name"))
all_genes.st <- as.character(res_tab.st$gene_id)

sig_num.st <- res_tab.st %>%
  filter(padj < 0.001) %>%
  pull(rowname) %>% length()
# 1082 genes are differentially expressed in PF tumors vs other PAs
top15_genes.st <- res_tab.st %>%
  filter(log2FoldChange > 0) %>%
  filter(rowname %in% protein_genes) %>%
  arrange(padj) %>%
  top_n(n = -15, padj) %>% pull(rowname)

sigOE_st <- dplyr::filter(res_tab.st, padj < 0.05 & log2FoldChange > 0) %>%
  dplyr::select(gene_id) %>% na.omit()
sigOE_st <- as.character(sigOE_st$gene_id)

sigUE_st <- dplyr::filter(res_tab.st, padj < 0.05 & log2FoldChange < 0) %>%
  dplyr::select(gene_id) %>% na.omit()
sigUE_st <- as.character(sigUE_st$gene_id)

stUE_GO <- enrichGO(gene = sigUE_st, 
                    universe = all_genes.st,
                    OrgDb = org.Hs.eg.db, 
                    keyType = 'ENSEMBL',
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE) %>%
  enrichplot::pairwise_termsim() %>%
  dotplot(x = "GeneRatio",
          showCategory = 20,
          font.size = 8,
          title = "Pathways under-enriched in supratentorial tumors")
stOE_GO <- enrichGO(gene = sigOE_st, 
                    universe = all_genes.st,
                    OrgDb = org.Hs.eg.db, 
                    keyType = 'ENSEMBL',
                    ont = "BP", 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = 0.05, 
                    readable = TRUE) %>%
  enrichplot::pairwise_termsim() %>%
  dotplot(x = "GeneRatio",
          showCategory = 20,
          font.size = 8,
          title = "Pathways enriched in supratentorial tumors")

pdf(paste0(fwd, "st_GO.pdf"))
par(mfrow = c(1, 2))
print(stUE_GO)
print(stOE_GO)
dev.off()

st_gsea_file <- res_tab.st %>%
  filter(padj < 0.01) %>%
  select(rowname, log2FoldChange) %>%
  arrange(-log2FoldChange)
write.table(x = st_gsea_file, 
            file = paste0(wd, "PA/PA_Data/st_geneList.txt"),
            sep = "\t",
            row.names = FALSE)

st_gProf_list <- res_tab.st %>%
  filter(padj < 0.05 & log2FoldChange >= 0) %>%
  pull(rowname)
st_gProf_list[st_gProf_list %in% protein_genes] -> st_gProf_list
write(x = st_gProf_list,
      file = paste0(wd, "PA/PA_Data/st_OE_geneList.txt"),
      sep = "\t")

save(annotations_ahb, dds_pca, dds_sp, dds_pf, dds_st, 
      file = paste0(wd, "PA/PA_Data/PiloA_GSEA_13.02.2022.RData"))


library(edgeR)
logCPM <- cpm(count_data, prior.count = 2, log = TRUE)
zScore <- t(scale(t(logCPM)))
genes_for_fig <- c(top15_genes.pf[1:10], top15_genes.spinal[1:10], top15_genes.st[1:10])
genes_for_fig <- genes_for_fig[genes_for_fig %in% rownames(dataScale)]
mat_for_ht <- zScore[genes_for_fig,]

sample_order_by_loc <- covariate_data %>%
  arrange(Location) %>%
  rownames()
mat_for_ht[,sample_order_by_loc] -> mat_for_ht
col_fun = colorRamp2(c(-1.5, 0, 1.5), c("blue", "white", "red"))
ht_list <- Heatmap(mat_for_ht,
                   border = "black",
                   col = col_fun,
                   show_column_names = FALSE,
                   cluster_rows = FALSE,
                   cluster_columns = FALSE,
                   heatmap_legend_param = list(
                     title = "Relative expression",
                     direction = "horizontal",
                     at = c(-2.5, 0, 2.5),
                     border = "black",
                     legend_width = unit(6, "cm"),
                     title_position = "topcenter"
                   ))
draw(ht_list, heatmap_legend_side = "bottom")

########## Single sample GSEA ##########
source(file = "prinses_max/PiloA/ssGSEA.R")
library(readxl)
library(DESeq2)

fid <- "norm_counts.gct" 
writeLines(c("#1.2", paste(nrow(counts(dds_pca, normalized=TRUE)), ncol(counts(dds_pca, normalized=TRUE)) - 2, collapse="\t")), fid, sep="\n")
write.table(counts(dds_pca, normalized=TRUE), file=fid, quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t", append = TRUE)

genesets <- read_excel("PA/Analysis/GSEA/Immune_metagenes.xlsx")
genesets <- genesets %>%
  select(Metagene, `Cell type`)
celltypes <- unique(genesets$`Cell type`)
celltype_genes <- list()
for (i in 1:length(celltypes)) {
  immune_cell <- celltypes[i]
  celltype_genes[[immune_cell]] <- genesets$Metagene[genesets$`Cell type` == immune_cell]
}

logcounts <- log2(counts(dds_pca, normalized=TRUE) + 1)
system.time(assign('ssGSEA_res', ssgsea(logcounts, celltype_genes, scale = FALSE, norm = FALSE)))
ssGSEA_res.scale <- (ssGSEA_res - rowMeans(ssGSEA_res))/(rowSds(as.matrix(ssGSEA_res)))[row(ssGSEA_res)]
dist <- cor(ssGSEA_res, method = "kendall")

ssGSEA_ht <- Heatmap(ssGSEA_res.scale,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        top_annotation = HeatmapAnnotation(Location = covariate_data$Location,
                                           col = col_annotations,
                                           show_legend = TRUE),
        show_row_names = TRUE,
        show_column_names = FALSE,
        show_row_dend = TRUE,
        heatmap_legend_param = list(
          title = "Relative expression",
          direction = "horizontal",
          at = c(-3, 0, 3),
          border = "black",
          legend_width = unit(3, "cm"),
          title_position = "topcenter"))

pdf(paste0(fwd, "ssGSEA_clusters.pdf"), width = 9, height = 8)
draw(ssGSEA_ht, heatmap_legend_side = "bottom")
dev.off()

########## Immune Panel ##########
immune_genes <- read_excel("PA/PA_Data/LBL-10043-08_nCounter_PanCancer_Immune_Profiling_Panel_Gene_List.xlsx", 
                                                                               sheet = "Annotations", skip = 1)
immune_genes <- immune_genes$`Gene Name`[1:770]

countsLog_immune <- log(cpm_data[immune_genes,] + 1)
varGenes_immune <- apply(countsLog_immune, 1, var)
sig_immuneGenes <- names(varGenes_immune[which(varGenes_immune > 0.75)])

mat <- dataScale[sig_immuneGenes[sig_immuneGenes %in% rownames(dataScale)],] %>% as.matrix()

immuneMarker_ht <- Heatmap(mat,
                           clustering_distance_columns = "kendall",
                           clustering_distance_rows = "kendall",
                           top_annotation = HeatmapAnnotation(Location = covariate_data$Location,
                                                              col = col_annotations, 
                                                              show_legend = TRUE),
                           show_row_names = TRUE,
                           show_column_names = FALSE,
                           show_row_dend = TRUE,
                           heatmap_legend_param = list(
                             title = "Relative expression",
                             direction = "horizontal",
                             at = c(-3, 0, 3),
                             border = "black",
                             legend_width = unit(3, "cm"),
                             title_position = "topcenter"))

pdf(paste0(fwd, "immuneMarker_clusters.pdf"), width = 9, height = 8)
draw(immuneMarker_ht, heatmap_legend_side = "bottom")
dev.off()


########## Deconvolution analysis ##########
### Make bulk and single-cell expression sets
scPA_eset <- readRDS("PA/Analysis/Deconvolution/DWLS/PA_singleCellReitman.RDS")
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
  mutate(Sample = as.character(rownames(deconv_res))) 
%>%
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

########## CIBERSORTx ##########
# write.csv(count_data, file = paste0(wd, "PA/Analysis/Deconvolution/CIBERSORTx/Bulk_PA_forDeconv.csv"))

PA_cibersortx_LM22 <- read.csv("PA/Analysis/Deconvolution/CIBERSORTx/PA_cibersortx_LM22.csv")
PA_lm22 <- PA_cibersortx_LM22[,1:23] %>%
  column_to_rownames(var = "Mixture") %>%
  dplyr::select(!Dendritic.cells.resting) %>%
  as.matrix() %>% t()
covariate_lm22 <- covariate_data[intersect(rownames(covariate_data), colnames(PA_lm22)),]
PA_lm22 <- PA_lm22[,intersect(rownames(covariate_data), colnames(PA_lm22))]


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

lm22_ht <- Heatmap(PA_lm22,
                           clustering_distance_columns = "kendall",
                           clustering_distance_rows = "manhattan",
                           top_annotation = HeatmapAnnotation(Location = covariate_data$Location,
                                                              col = col_annotations, 
                                                              show_legend = TRUE),
                           show_row_names = TRUE,
                           show_column_names = FALSE,
                           show_row_dend = TRUE,
                           heatmap_legend_param = list(
                             title = "Cell type composition (%)",
                             direction = "horizontal",
                             at = c(0, 1.5),
                             border = "black",
                             legend_width = unit(3, "cm"),
                             title_position = "topcenter"))

pdf(paste0(fwd, "cibersort_clusters.pdf"), width = 9, height = 8)
draw(lm22_ht, heatmap_legend_side = "bottom")
dev.off()
