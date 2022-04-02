# Pilocytic astrocytoma transcriptomic analysis
# Last updated: 1/04/2022

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
               "Posterior Fossa" = "#636EFA", 
               "Supratentorial" = "#00CC96"),
  LocationSpecific = c("cerebellum/4th" = "#636EFA",
                       "spinal cord" = "#EF553B",
                       "brain stem" = "black",
                       "diencephalon" = "#00CC96",
                       "cerebral hemisphere" = "yellow"),
  Age = colorRamp2(c(0,20), c("#edf8e9", "#74c476"))
)

wd <- "/Users/jrozowsky/Library/Mobile Documents/com~apple~CloudDocs/Documents/PMC/"
setwd(wd)
fwd <- paste0(wd, "PA/Analysis/Figures/") # figure working directory

annotations_ahb <- read.csv("PA/Analysis/GSEA/annotations_ahb.csv")
annotations_ahb %>% filter(gene_biotype == "protein_coding") %>% pull(gene_name) -> protein_genes

########## Prepare Data ##########
# compile all patient IDs
file.path <- paste0(wd, "Data/PMC/RNA_jacob/")
files <- paste0(list.files(path = file.path))

meta_data <- read_excel(path = "PA/Cohort/ClinicalData_updates.xlsx", sheet = "CohortList")

patient_counts <- data.frame()
patient_id_list <- c()
for (i in 1:length(files)) {
  print(i)
  sample_id <- files[i] %>% strsplit(split = "_") %>% sapply(getElement, 1)
  print(sample_id)
  if (meta_data$UniqueSample[meta_data$PMABM == sample_id] == "Yes") {
    patient_id <- meta_data %>% 
      filter(UniqueSample == "Yes") %>%
      filter(PMABM == sample_id) %>%
      pull(`Subject ID`)
    print(patient_id)
    
    if (!patient_id %in% patient_id_list) {
      patient_id_list <- c(patient_id_list, patient_id)
      patient_file <- paste0(file.path, files[i])
      colnames <- c("ID", "Counts", "CPM", "FPKM", "Chr", "Start", "End", "Strand", "Length", "GeneID", "GeneName", "TranscriptID")
      patient_counts <- read.table(file = patient_file, header = FALSE)
      colnames(patient_counts) <- colnames
      if (patient_counts$ID[1] == "ID") {
        patient_counts <- patient_counts[-1,]
      }
      patient_counts <- patient_counts %>%
        select(GeneName, Counts)
      patient_counts$Counts <- as.numeric(patient_counts$Counts)
      # there are multiple rows to each gene name. add up counts so that there is one row per gene
      patient_counts <- aggregate(Counts ~ GeneName, patient_counts, sum)
      
      colnames(patient_counts)[2] <- patient_id
      
      if (i == 1) {
        genes <- patient_counts$GeneName
        pa_counts <- data.frame(patient_counts)
        rm(patient_counts)
      } else {
        pa_counts <- merge(pa_counts, patient_counts, by = "GeneName")
        rm(patient_counts)
      }
    }
  }
}

########## Analyze Data ##########
load("PA/Cohort/Cohort_BulkRNAseq.RData")

pa_counts <- pa_counts %>% column_to_rownames(var = "GeneName")
covariate_data <- meta_data %>% 
  filter(UniqueSample == "Yes") %>%
  column_to_rownames(var = "Subject ID") %>%
  select("Location", "LocationSpecific", "AgeMonths", "Sex", "Molecular.1", "Molecular.2")
covariate_data$Location <- as.factor(covariate_data$Location)
covariate_data$LocationSpecific <- as.factor(covariate_data$LocationSpecific)
covariate_data$Sex <- as.factor(covariate_data$Sex)
unique_samples <- rownames(covariate_data)

xy_genes <- read.delim2("PA/PA_Data/xy_genes.txt", header = TRUE, sep = "\t", dec = ".")

## remove x and y linked genes
genes_to_use <- setdiff(rownames(pa_counts), xy_genes$Approved.symbol)
pa_cpm <- apply(pa_counts, 2, function(x) (x/sum(x))*1000000)
countsLog <- log(pa_cpm[genes_to_use,] + 1)
varGenes <- apply(countsLog, 1, var)
meanGenes <- apply(countsLog, 1, mean)
nFeatures = 5000
varFeatures <- names(varGenes)[order(varGenes, decreasing = T)][c(1:nFeatures)]
dataScale <- apply(countsLog, 2, function(x) (x - meanGenes)/varGenes)

########## Whole transcriptome analysis ##########
load(file = "PA/Cohort/CohortBulkRNAseq_01.04.2022.RData")
### Dimensionality reduction ###
### PCA ###
library(DESeq2)
order <- colnames(pa_counts)
covariate_data <- covariate_data[order,]
# change medulla to brainstem
covariate_data$LocationSpecific <- gsub(x = covariate_data$LocationSpecific, pattern = "medulla", replacement = "brain stem")
dds <- DESeqDataSetFromMatrix(countData = pa_counts,
                              colData = covariate_data,
                              design = ~ 1)
dds$LocationSpecific <- factor(dds$LocationSpecific)
dds$Location <- factor(dds$Location)
design(dds) <- ~LocationSpecific
dds <- DESeq(dds, test = "LRT", reduced = ~ 1)
vsd <- vst(dds, blind = FALSE)

save(pa_counts,
     meta_data,
     dataScale,
     varFeatures,
     covariate_data,
     dds,
     file = paste0(wd, "PA/PA_Data/CohortBulkRNAseq_02.04.2022.RData"))
# Select the most variably expressed genes -> these are non-x/y genes
### UMAP ###
library(umap)
library(ggforce)
numGenes <- 600
seed <- 1000

pca <- prcomp(t(assay(vsd)[varFeatures[1:numGenes], ]))
set.seed(seed); umap_res <- umap(pca$x, alpha = 0.62, gamma = 1)
umap_data <- cbind(umap_res$layout, data.frame(dds_pca@colData))
bulkUmapPlot <- ggplot(data = umap_data, aes(x = `1`, y = `2`)) +
  geom_point(aes(fill = LocationSpecific), colour = "grey30", pch = 21, size = 2) + 
  xlab("UMAP1") + ylab("UMAP2") +
  guides(col = guide_legend(ncol = 1)) +
  scale_fill_manual(name = "Location", values = col_annotations$LocationSpecific) +
  theme_classic() +
  labs(subtitle = paste0("UMAP: varFeatures = ", numGenes, "; set.seed = ", seed)) +
  theme(panel.background = element_rect(colour = "grey30", size=1),
        axis.line = element_blank(),
        plot.title = element_text(size = 20),
        axis.ticks = element_blank(),
        axis.text = element_blank())

pdf(paste0(fwd, "BulkCohortUMA_01.04.2022.pdf"), width = 6.5, height = 5)
print(bulkUmapPlot)
dev.off()

### Unsupervised clustering ###
library(ComplexHeatmap)
# Identify protein-coding genes
varFeatures_coding <- varFeatures[varFeatures %in% protein_genes]

# calculate inter-patient correlation based on protein-coding gene expression
# corr_res <- cor(dataScale[varFeatures_coding[1:262],], method = "pearson")
Heatmap(assay(vsd)[varFeatures[1:numGenes], ],
        show_heatmap_legend = TRUE,
        clustering_distance_rows = "euclidean",
        clustering_distance_columns = "euclidean",
        top_annotation = HeatmapAnnotation(Sex = covariate_data$Sex,
                                           LocationSpecific = covariate_data$LocationSpecific,
                                           col = col_annotations, 
                                           show_legend = TRUE),
        show_row_names = FALSE,
        show_column_names = FALSE,
        show_row_dend = FALSE)
pdf(paste0(fwd, "heatmap.pdf"))
print(heatmap)
dev.off()

########## Functional enrichment analysis ##########
library(fgsea)
library(org.Hs.eg.db)
library(fgsea)
library(clusterProfiler)
library(DESeq2)
load("PA/PA_Data/CohortBulkRNAseq_02.04.2022.RData")

spinal_patients <- rownames(covariate_data)[covariate_data$LocationSpecific == "spinal cord"]
posteriorFossa_patients <- rownames(covariate_data)[covariate_data$LocationSpecific %in% c("cerebellum/4th", "brain stem")]
diencephalon_patients <- rownames(covariate_data)[covariate_data$LocationSpecific == "diencephalon"]
cerebralHemisphere_patients <- rownames(covariate_data)[covariate_data$LocationSpecific == "cerebral hemisphere"]

covariate_data <- covariate_data %>%
  mutate(Spinal = as.factor(ifelse(rownames(covariate_data) %in% spinal_patients, "Spinal", "not.Spinal")),
         PosteriorFossa = as.factor(ifelse(rownames(covariate_data) %in% posteriorFossa_patients, "PosteriorFossa", "not.PosteriorFossa")),
         Diencephalon = as.factor(ifelse(rownames(covariate_data) %in% diencephalon_patients, "Diencephalon", "not.Diencephalon")),
         Cerebral = as.factor(ifelse(rownames(covariate_data) %in% cerebralHemisphere_patients, "CerebralHemisphere", "not.CerebralHemisphere")))

# Make all DDS objects
dds_sp <- DESeqDataSetFromMatrix(countData = pa_counts,
                                    colData = covariate_data,
                                    design = ~ Spinal)
dds_sp <- DESeq(dds_sp)

dds_pf <- DESeqDataSetFromMatrix(countData = pa_counts,
                                 colData = covariate_data,
                                 design = ~ PosteriorFossa)
dds_pf <- DESeq(dds_pf)

dds_d <- DESeqDataSetFromMatrix(countData = pa_counts,
                                 colData = covariate_data,
                                 design = ~ Diencephalon)
dds_d <- DESeq(dds_d)

dds_ch <- DESeqDataSetFromMatrix(countData = pa_counts,
                                 colData = covariate_data,
                                 design = ~ Cerebral)
dds_ch <- DESeq(dds_ch)

save(dds_ch, dds_d, dds_pf, dds_sp, file = paste0(wd, "PA/PA_Data/DESeqBulkCohort_02.04.2022.RData"))

cnsLocation <- c("Cerebral", "Diencephalon", "PosteriorFossa", "Spinal")
contrasts <- list(Cerebral = c("Cerebral", "CerebralHemisphere", "not.CerebralHemisphere"),
                  Diencephalon = c("Diencephalon", "Diencephalon", "not.Diencephalon"),
                  PosteriorFossa = c("PosteriorFossa", "PosteriorFossa", "not.PosteriorFossa"),
                  Spinal = c("Spinal", "Spinal", "not.Spinal"))
dds_objects <- objects()[grep("dds_", objects())]
for (i in 1:length(dds_objects)) {
  # initialize the location and datasets
  dds_loc <- get(dds_objects[i])
  temp <- as.character(dds_loc@design) %>% strsplit(split = "~") %>% sapply(getElement, 1)
  location <- temp[2]
  contrast <- contrasts[[location]]
  print(location)
  print(contrast)
  
  # perform differential expression
  res <- results(dds_loc, alpha = 0.05, test = "Wald", contrast = contrast)
  res.tab <- data.frame(res[order(res$padj)[1:5000],])
  res.tab <- left_join(res.tab %>% rownames_to_column(), annotations_ahb, by=c("rowname"="gene_name"))
  
  # select most enriched genes for this location
  nGenes <- 15
  top_genes <- res.tab %>%
    filter(rowname %in% protein_genes) %>%
    top_n(n = nGenes, log2FoldChange) %>%
    arrange(-log2FoldChange) %>% pull(rowname)
  
  OEgenes <- filter(res.tab, padj < 0.05 & log2FoldChange > 0) %>% arrange(-log2FoldChange) %>% pull(gene_id)
  UEgenes <- filter(res.tab, padj < 0.05 & log2FoldChange < 0) %>% arrange(log2FoldChange) %>% pull(gene_id)
  
  UE_GO <- enrichGO(gene = UEgenes,
                      universe = as.character(res.tab$gene_id),
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
            title = paste0("Pathways under-enriched in ", location, " tumors"))
  OE_GO <- enrichGO(gene = OEgenes,
           universe = as.character(res.tab$gene_id),
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
            title = paste0("Pathways over-enriched in ", location, " tumors"))

  pdf(paste0(fwd, "Bulk.GOenrichmentPlot.", location, "_02.04.2022.pdf"))
  par(mfrow = c(1, 2))
  print(UE_GO)
  print(OE_GO)
  dev.off()
  
  gProf.genes_entrz <- filter(res.tab, padj < 0.05 & log2FoldChange > 0) %>% filter(gene_biotype == "protein_coding") %>% arrange(-log2FoldChange) %>% pull(gene_id)
  gProf.genes_symbol <- filter(res.tab, padj < 0.05 & log2FoldChange > 0) %>% filter(gene_biotype == "protein_coding") %>% arrange(-log2FoldChange) %>% pull(rowname)
  write(x = gProf.genes_entrz,
        file = paste0(wd, "PA/PA_Data/gProfiler_entrz.", location, "_02.04.2022.txt"),
        sep = "\t")
  write(x = gProf.genes_symbol,
        file = paste0(wd, "PA/PA_Data/gProfiler_symbol.", location, "_02.04.2022.txt"),
        sep = "\t")
}

library(edgeR)
logCPM <- cpm(pa_counts, prior.count = 2, log = TRUE)
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
meta_genes <- unique(intersect(sc_genes, rownames(pa_counts)))

bulkPA_eset <- ExpressionSet(assayData = as.matrix(pa_counts[meta_genes,]))

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
# write.csv(pa_counts, file = paste0(wd, "PA/Analysis/Deconvolution/CIBERSORTx/Bulk_PA_forDeconv.csv"))

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
