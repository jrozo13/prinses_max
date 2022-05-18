# Pilocytic astrocytoma transcriptomic analysis: Functional Enrichment
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
  Location = c("spinal" = "#EF553B", 
               "posterior fossa" = "#636EFA", 
               "supratentorial" = "#00CC96"),
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

########## Make DEG Objects ##########
library(fgsea)
library(DESeq2)
load("PA/PA_Data/CohortBulkRNAseq_04.04.2022.RData")

levels(covariate_data$Location) <- c("spinal", "posterior fossa", "supratentorial")
posteriorFossa_patients <- rownames(covariate_data)[covariate_data$Location == "posterior fossa"]
supratentorial_patient <- rownames(covariate_data)[covariate_data$Location == "supratentorial"]
spinal_patients <- rownames(covariate_data)[covariate_data$Location == "spinal"]

covariate_data <- covariate_data %>%
  mutate(PosteriorFossa = as.factor(ifelse(rownames(covariate_data) %in% posteriorFossa_patients, "PosteriorFossa", "not.PosteriorFossa")),
         Supratentorial = as.factor(ifelse(rownames(covariate_data) %in% supratentorial_patient, "Supratentorial", "not.Supratentorial")),
         Spinal = as.factor(ifelse(rownames(covariate_data) %in% spinal_patients, "Spinal", "not.Spinal")))

dds_pf <- DESeqDataSetFromMatrix(countData = pa_counts,
                                 colData = covariate_data,
                                 design = ~ PosteriorFossa)
dds_pf <- DESeq(dds_pf)

dds_st <- DESeqDataSetFromMatrix(countData = pa_counts,
                                colData = covariate_data,
                                design = ~ Supratentorial)
dds_st <- DESeq(dds_st)

dds_sp <- DESeqDataSetFromMatrix(countData = pa_counts,
                                 colData = covariate_data,
                                 design = ~ Spinal)
dds_sp <- DESeq(dds_sp)

#save(dds_pf, dds_st, dds_sp, file = paste0(wd, "PA/PA_Data/DESeqBulkCohort.3Location_12.05.2022.RData"))

########## Differentially expressed genes ##########
load("PA/PA_Data/DESeqBulkCohort.3Location_12.05.2022.RData")
load("PA/PA_Data/CohortBulkRNAseq_04.04.2022.RData")
library(DESeq2)
library(org.Hs.eg.db)
library(fgsea)
library(clusterProfiler)

cnsLocation <- c("PosteriorFossa", "Supratentorial", "Spinal")
contrasts <- list(PosteriorFossa = c("PosteriorFossa", "PosteriorFossa", "not.PosteriorFossa"),
                  Supratentorial = c("Supratentorial", "Supratentorial", "not.Supratentorial"),
                  Spinal = c("Spinal", "Spinal", "not.Spinal"))
dds_objects <- c("dds_pf", "dds_st", "dds_sp")

for (i in 1:length(dds_objects)) {
  dds_loc <- get(dds_objects[i])
  temp <- as.character(dds_loc@design) %>% strsplit(split = "~") %>% sapply(getElement, 1)
  location <- temp[2]
  contrast <- contrasts[[location]]
  print(location)
  print(contrast)
  
  # perform differential expression
  res <- results(dds_loc, alpha = 0.05, test = "Wald", contrast = contrast)
  res.tab <- data.frame(res[order(res$padj),][1:5000,])
  res.tab <- left_join(res.tab %>% rownames_to_column(), annotations_ahb, by=c("rowname"="gene_name"))
  
  # select most enriched genes for this location
  nGenes <- 15
  top_genes <- res.tab %>%
    filter(rowname %in% protein_genes) %>%
    filter(log2FoldChange > 0) %>%
    arrange(padj) %>%
    pull(rowname) %>%
    .[1:nGenes]
  assign(paste0(location, "_genes"), top_genes)

  OEgenes <- filter(res.tab, padj < 0.05 & log2FoldChange > 0) %>% arrange(-log2FoldChange) %>% pull(gene_id) %>% unique()
  UEgenes <- filter(res.tab, padj < 0.05 & log2FoldChange < 0) %>% arrange(log2FoldChange) %>% pull(gene_id)%>% unique()
  
  # UE_GO <- enrichGO(gene = UEgenes,
  #                   universe = as.character(res.tab$gene_id),
  #                   OrgDb = org.Hs.eg.db,
  #                   keyType = 'ENSEMBL',
  #                   ont = "BP",
  #                   pAdjustMethod = "BH",
  #                   qvalueCutoff = 0.05,
  #                   readable = TRUE) %>%
  #   enrichplot::pairwise_termsim() %>%
  #   dotplot(x = "GeneRatio",
  #           showCategory = 20,
  #           font.size = 6) + 
  #   ggtitle(label = paste0("Pathways under-enriched in ", location, " tumors")) +
  #   theme(plot.title = element_text(size = 10))
  # OE_GO <- enrichGO(gene = OEgenes,
  #                   universe = as.character(res.tab$gene_id),
  #                   OrgDb = org.Hs.eg.db,
  #                   keyType = "ENSEMBL",
  #                   ont = "BP",
  #                   pAdjustMethod = "BH",
  #                   qvalueCutoff = 0.05,
  #                   readable = TRUE) %>%
  #   enrichplot::pairwise_termsim() %>%
  #   dotplot(x = "GeneRatio",
  #           showCategory = 20,
  #           font.size = 6) + 
  #   ggtitle(label = paste0("Pathways over-enriched in ", location, " tumors")) +
  #   theme(plot.title = element_text(size = 10))

  # pdf(paste0(fwd, "Bulk.GOenrichmentPlot.", location, "_12.05.2022.pdf"), width = 10, height = 5)
  # grid.arrange(UE_GO, OE_GO, ncol = 2)
  # dev.off()
  
  gProf.genes_entrz <- filter(res.tab, padj < 0.05 & log2FoldChange > 0) %>% filter(gene_biotype == "protein_coding") %>% arrange(-log2FoldChange) %>% pull(gene_id)
  gProf.genes_symbol <- filter(res.tab, padj < 0.05 & log2FoldChange > 0) %>% filter(gene_biotype == "protein_coding") %>% arrange(-log2FoldChange) %>% pull(rowname)
  # write(x = gProf.genes_entrz,
  #       file = paste0(wd, "PA/PA_Data/gProfiler_entrz.", location, "_12.05.2022.txt"),
  #       sep = "\t")
  # write(x = gProf.genes_symbol,
  #       file = paste0(wd, "PA/PA_Data/gProfiler_symbol.", location, "_12.05.2022.txt"),
  #       sep = "\t")
}

de_genes <- c(PosteriorFossa_genes, Supratentorial_genes, Spinal_genes)
label_genes <- c(c("GRM4", "DAO", "SPX", "CDH15", "EN2", "PLEKHD1"),
                 c("FOXG1", "CCKAR", "SIX3", "LRRC7", "SLIT1"),
                 c("HOXB7", "HOXB5", "HOXC8", "HOXC10", "HOXC5"))
label_index <- match(label_genes, de_genes)
sample_order <- c(posteriorFossa_patients, supratentorial_patient, spinal_patients)
annotation_order <- covariate_data[sample_order,]
annotation_order$Location <- factor(annotation_order$Location, levels = c("posterior fossa", "supratentorial", "spinal"))

library(edgeR)
library(grid)
library(ComplexHeatmap)
library(gridtext)
pa_cpm <- cpm(pa_counts, prior.count = 2, log = TRUE)
pa_zScore <- t(scale(t(pa_cpm)))

pa_zScore <- pa_zScore[de_genes, sample_order]

col_fun = colorRamp2(c(-0.75, 0, 1), c("navy", "black", "red"))
# col_fun = colorRamp2(c(-1, 0, 1.5), c("navy", "white", "red"))
geneByLoc_Figure <- grid.grabExpr(draw(Heatmap(
  matrix = pa_zScore,
  col = col_fun,
  
  cluster_rows = FALSE,
  row_split = c(rep(1,nGenes), rep(2,nGenes), rep(3,nGenes)),
  row_title = NULL,
  show_row_names = FALSE,
  
  cluster_columns = FALSE, 
  column_split = annotation_order$Location,
  show_column_names = FALSE, 
  column_title = NULL,
  show_heatmap_legend = FALSE,
  right_annotation = rowAnnotation(link = anno_mark(at = label_index,
                                                    labels = label_genes,
                                                    labels_gp = gpar(fontface = "italic",
                                                                     fontsize=15))),
  
  # top_annotation = HeatmapAnnotation(Location = annotation_order$Location,
  #                                    col = col_annotations, 
  #                                    show_legend = FALSE,
  #                                    show_annotation_name = FALSE)
  )),
  width = 7, height = 7)

geneByLoc_Legend <- grid.grabExpr(draw(Legend(col_fun = col_fun, 
                                              title = "Relative Expression",
                                              direction = "horizontal",
                                              at = c(-2, 0, 2),
                                              legend_width = unit(4, "cm"),
                                              title_position = "topcenter")),
                                  width = 3)

pdf(paste0(fwd, "Bulk.DEGs_Location.pdf"), width = 9, height = 5)
plot_grid(geneByLoc_Figure)
dev.off()

pdf(paste0(fwd, "Bulk.DEGs_Location.Legend.pdf"), width = 3, height = 3)
plot_grid(geneByLoc_Legend)
dev.off()

########## Functional Enrichment ##########
load("PA/PA_Data/PiloA_Data_09.02.2022.RData")
load("PA/PA_Data/PiloA_GSEA_13.02.2022.RData", oldData <- new.env())
dds_pf.old <- oldData$dds_pf
dds_st.old <- oldData$dds_st
dds_sp.old <- oldData$dds_sp

dds_objects.old <- c("dds_pf.old", "dds_st.old", "dds_sp.old")
cnsLocation.old <- c("PosteriorFossa", "Supratentorial", "Spinal")
contrasts.old <- list(PF_Loc = c("PF_Loc", "Posterior fossa", "Not.PF"),
                      ST_Loc = c("ST_Loc", "Supratentorial", "Not.ST"),
                      Spinal_Loc = c("Spinal_Loc", "Spinal", "Not.spinal"))
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
for (i in 1:length(dds_objects.old)) {
  dds_loc <- get(dds_objects[i])
  temp <- as.character(dds_loc@design) %>% strsplit(split = "~") %>% sapply(getElement, 1)
  location <- temp[2]
  contrast <- contrasts.old[[location]]
  print(location)
  print(contrast)
  
  # perform differential expression
  res <- results(dds_loc, alpha = 0.05, test = "Wald", contrast = contrast)
  res.tab <- data.frame(res[order(res$padj),][1:5000,])
  res.tab <- left_join(res.tab %>% rownames_to_column(), annotations_ahb, by=c("rowname"="gene_name"))
  
  OEgenes <- filter(res.tab, padj < 0.05 & log2FoldChange > 0) %>% arrange(-log2FoldChange) %>% pull(gene_id) %>% unique()
  UEgenes <- filter(res.tab, padj < 0.05 & log2FoldChange < 0) %>% arrange(log2FoldChange) %>% pull(gene_id)%>% unique()
  
  UE_GO <- enrichGO(gene = UEgenes,
                    universe = as.character(res.tab$gene_id),
                    OrgDb = org.Hs.eg.db,
                    keyType = 'ENSEMBL',
                    ont = "BP",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE) %>%
    enrichplot::pairwise_termsim()
  UE_GO <- UE_GO@result %>%
    filter(p.adjust < 0.05)
  UE_GO$Location_Enrich <- paste0(location, ".down")
  assign(paste0(location, ".down"), UE_GO)

  OE_GO <- enrichGO(gene = OEgenes,
                    universe = as.character(res.tab$gene_id),
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENSEMBL",
                    ont = "BP",
                    pAdjustMethod = "BH",
                    qvalueCutoff = 0.05,
                    readable = TRUE) %>%
    enrichplot::pairwise_termsim()
  OE_GO <- OE_GO@result %>%
    filter(p.adjust < 0.05)
  OE_GO$Location_Enrich <- paste0(location, ".up")
  assign(paste0(location, ".up"), OE_GO)
}

pf.up.geneset <- c("neurotransmitter transport",
                   "nervous system development",
                   "cerebellum development",
                   "dendrite morphogenesis")
st.up.geneset <- c("tRNA metabolic process",
                   "adaptive immune response",
                   "regulation of T cell differentiation",
                   "mitochondrial translation")
spinal.up.geneset <- c("embryonic skeletal system development",
                       "RNA processing",
                       "chromosome organization",
                       "cranial nerve development")

pf.enrich <- PF_Loc.up %>% filter(Description %in% pf.up.geneset) %>% dplyr::select(Description, p.adjust) %>%
  ggplot(aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat = "identity", fill = col_annotations$Location[["posterior fossa"]]) +
  coord_flip() +
  geom_text(aes(label = Description,
                y = 4/30),
            hjust = "left") +
  ylab("-log10(p-adj)") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        legend.position = "none")

st.enrich <- ST_Loc.up %>% filter(Description %in% st.up.geneset) %>% dplyr::select(Description, p.adjust) %>%
  ggplot(aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat="identity", fill = col_annotations$Location[["supratentorial"]]) +
  coord_flip() +
  geom_text(aes(label = Description,
                y = 4/30),
            hjust = "left") +
  ylab("-log10(p-adj)") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        legend.position = "none")

sp.enrich <- Spinal_Loc.up %>% filter(Description %in% spinal.up.geneset) %>% dplyr::select(Description, p.adjust)
sp.enrich[1,1] <- "embryonic skeletal development"
sp.enrich <- sp.enrich %>%
  ggplot(aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust))) +
  geom_bar(stat="identity", fill = col_annotations$Location[["spinal"]]) +
  coord_flip() +
  geom_text(aes(label = Description,
                y = 7/30),
            hjust = "left") +
  ylab("-log10(p-adj)") +
  theme_classic() +
  theme(axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 8),
        axis.title.x = element_text(size = 10),
        legend.position = "none")


pdf(paste0(fwd, "Bulk.FunctionEnrichDEGs.Location.pdf"), width = 12, height = 7)
plot_grid(geneByLoc_Figure,
          plot_grid(pf.enrich, st.enrich, sp.enrich, nrow = 3),
          ncol = 2,
          rel_widths = c(10,3))
dev.off()

plot_grid(geneByLoc_Figure)

