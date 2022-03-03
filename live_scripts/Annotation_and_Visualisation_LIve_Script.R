library(AnnotationHub)
library(AnnotationDbi)
library(ensembldb)
library(DESeq2)
library(tidyverse)

# BiocManager::install("ensembldb")

ddsObj.interaction <- readRDS("RObjects/DESeqDataSet.interaction.rds")
results.interaction.d11 <- readRDS("RObjects/DESeqResults.interaction_d11.rds")
results.interaction.d33 <- readRDS("RObjects/DESeqResults.interaction_d33.rds")

# Start annotationhub

ah <- AnnotationHub()

query(ah, c("EnsDb", "Mus musculus", "102"))

MouseEnsDb <- query(ah, c("EnsDb", "Mus musculus", "102"))[[1]]

annotations <- genes(MouseEnsDb, return.type = "data.frame")
colnames(annotations)

annot <- annotations %>% 
  select(gene_id, gene_name, entrezid) %>% 
  filter(gene_id %in% rownames(results.interaction.d11))

length(unique(annot$entrezid))
sum(is.na(annot$entrezid))

# One we prepared ealier

ensemblAnnot <- readRDS("RObjects/Ensembl_annotations.rds")

annot.interaction.11 <- as.data.frame(results.interaction.d11) %>% 
  rownames_to_column("GeneID") %>% 
  left_join(ensemblAnnot, by = "GeneID") %>% 
  rename(logFC = log2FoldChange, FDR = padj)

write_tsv(annot.interaction.11, "results/Interaction.11_Results_Annotated.txt")

# Histogram p-values

hist(annot.interaction.11$pvalue)

# Shrinking the log2FoldChange

ddsShrink.11 <- lfcShrink(ddsObj.interaction,
                          res = results.interaction.d11,
                          type = "ashr")
shrinkTab.11 <- as.data.frame(ddsShrink.11) %>%
  rownames_to_column("GeneID") %>% 
  left_join(ensemblAnnot, by = "GeneID") %>% 
  rename(logFC = log2FoldChange, FDR = padj)

# MA plots

par(mfrow=c(1, 2))
plotMA(results.interaction.d11, alpha = 0.05)
plotMA(ddsShrink.11, alpha = 0.05)

# Volcano Plots

volcanoTab.11 <- shrinkTab.11 %>%
  mutate(`-log10(pvalue)` = -log10(pvalue))

ggplot(volcanoTab.11, aes(x = logFC, y = `-log10(pvalue)`)) +
  geom_point(aes(colour = FDR < 0.05), size = 1) +
  geom_text(data = ~top_n(.x, 1, wt=-FDR), aes(label = Symbol))

# Exercise 1 - Volcano plot for day 33

results.interaction.d33 <- readRDS("RObjects/DESeqResults.interaction_d33.rds")
ensemblAnnot <- readRDS("RObjects/Ensembl_annotations.rds")

# a) Shrink the log2 fold changes

ddsShrink.33 <- lfcShrink(ddsObj.interaction,
                          res = results.interaction.d33,
                          type = "ashr")

shrinkTab.33 <- as.data.frame(ddsShrink.33) %>% 
  rownames_to_column("GeneID") %>% 
  left_join(ensemblAnnot, by = "GeneID") %>% 
  rename(logFC = log2FoldChange, FDR = padj)

# b) Create -log10(pvalue) column

volcanoTab.33 <- shrinkTab.33 %>% 
  mutate(`-log10(pvalue)` = -log10(pvalue))

ggplot(volcanoTab.33, aes(x = logFC, y = `-log10(pvalue)`)) +
  geom_point(aes(colour = FDR < 0.05), size = 1)

# Exercise 2 - MA plot for day 33

maTab.33 <- shrinkTab.33 %>% 
  mutate(M = log2(baseMean))

ggplot(maTab.33, aes(x = M, y = logFC)) +
  geom_point(aes(colour = FDR < 0.05), size = 1) +
  scale_y_continuous(limit = c(-4, 4),
                     oob = scales::squish)

# Venn Diagram

library(ggvenn)

vennDat <- tibble(Geneid = rownames(results.interaction.d11)) %>% 
  mutate(Upregulated_11 = results.interaction.d11$padj < 0.05 &
                          !is.na(results.interaction.d11$padj) &
                          results.interaction.d11$log2FoldChange > 0) %>% 
  mutate(Downregulated_11 = results.interaction.d11$padj < 0.05 & 
           !is.na(results.interaction.d11$padj) & 
           results.interaction.d11$log2FoldChange < 0) %>%
  mutate(Upregulated_33 = results.interaction.d33$padj < 0.05 & 
           !is.na(results.interaction.d33$padj) & 
           results.interaction.d33$log2FoldChange > 0) %>%
  mutate(Downregulated_33 = results.interaction.d33$padj < 0.05 & 
           !is.na(results.interaction.d33$padj) & 
           results.interaction.d33$log2FoldChange < 0) 

ggvenn(vennDat, set_name_size = 4)

# Heatmap

library(ComplexHeatmap)
library(circlize)

sigGenes <- shrinkTab.11 %>%
  top_n(n = 300, wt = -FDR) %>% 
  pull("GeneID")

plotDat <- vst(ddsObj.interaction)[sigGenes,] %>% 
  assay()

z.mat <- t(scale(t(plotDat), center = TRUE, scale = TRUE))

# colour palette

myPalette <- c("royalblue3", "ivory", "orangered3")
myRamp <- colorRamp2(c(-2, 0, 2), myPalette)


Heatmap(z.mat, name = "z-score",
        col = myRamp,
        show_row_names = FALSE)

# create an annotation

ha1 <- HeatmapAnnotation(df = colData(ddsObj.interaction)[,c("Status", "TimePoint")])

Heatmap(z.mat, name = "z-score",
        col = myRamp,
        show_row_names = FALSE,
        split = 3,
        top_annotation = ha1,
        rect_gp = gpar(col = "lightgrey", lwd = 0.3))

ha1 <- HeatmapAnnotation(df = colData(ddsObj.interaction)[,c("Status", "TimePoint")],
                         col = list(Status = c("Uninfected" = "darkgreen",
                                               "Infected" = "palegreen"),
                                    TimePoint = c("d11" = "lightblue",
                                                  "d33" = "darkblue")))

Heatmap(z.mat, name = "z-score",
        col = myRamp,
        show_row_names = FALSE,
        split = 3,
        top_annotation = ha1,
        rect_gp = gpar(col = "lightgrey", lwd = 0.3))














