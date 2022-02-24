library(DESeq2)
library(tidyverse)

# reading in the data

txi <- readRDS("RObjects/txi.rds")
sampleinfo <- read_tsv("data/samplesheet_corrected.tsv", col_types = "cccc")

all(colnames(txi$counts)==sampleinfo$SampleName)

# making the deseq2 object

simple.model <- as.formula(~ Status)
model.matrix(simple.model, data = sampleinfo)

sampleinfo <- mutate(sampleinfo, Status = fct_relevel(Status, "Uninfected"))
model.matrix(simple.model, data = sampleinfo)

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = simple.model)

keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep,]

# DESeq2

# size factors

ddsObj <- estimateSizeFactors(ddsObj.filt)

normalizationFactors(ddsObj.filt)
normalizationFactors(ddsObj)

logcounts <- log2(counts(ddsObj, normalized = FALSE) + 1)

limma::plotMA(logcounts, array = 5, ylim = c(-5, 5))
abline(h=0, col="red")

logNormcounts <- log2(counts(ddsObj, normalized = TRUE) + 1)

limma::plotMA(logNormcounts, array = 5, ylim = c(-5, 5))
abline(h=0, col="red")

# estimate dispersion

ddsObj <- estimateDispersions(ddsObj)

plotDispEsts(ddsObj)

# wald test

ddsObj <- nbinomWaldTest(ddsObj)

# DESeq wrapper

ddsObj <- DESeq(ddsObj.filt)

# get results out

results.simple <- results(ddsObj, alpha = 0.05)
results.simple

# Exercise 1

sum(results.simple$padj < 0.05) 

sum(is.na(results.simple$padj))

sum(results.simple$padj < 0.05 & results.simple$log2FoldChange > 0, na.rm = TRUE) 

sum(results.simple$padj < 0.05 & results.simple$log2FoldChange < 0, na.rm = TRUE) 





