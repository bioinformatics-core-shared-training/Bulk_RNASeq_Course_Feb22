library(tximport)
library(DESeq2)
library(tidyverse)

# reading in sample data

sampleinfo <- read_tsv("data/samplesheet.tsv", col_types = c("cccc"))
sampleinfo %>%
  arrange(Status, TimePoint, Replicate)

# reading in the count data

files <- str_c("salmon/", sampleinfo$SampleName, "/quant.sf")
files <- set_names(files, sampleinfo$SampleName)

tx2gene <- read_tsv("references/tx2gene.tsv")

txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
str(txi)
head(txi$counts)

#Exercise 1

tpm <- tximport(files, 
                type = "salmon", 
                tx2gene = tx2gene, 
                countsFromAbundance = "lengthScaledTPM")

# make raw counts object

rawCounts <- round(txi$counts, 0)

# filtering

dim(rawCounts)

keep <- rowSums(rawCounts) > 5
table(keep, useNA = "always")

filtCounts <- rawCounts[keep,]
dim(filtCounts)

# transformations

summary(filtCounts)
boxplot(filtCounts, main = "Raw Counts", las = 2)

plot(rowMeans(filtCounts), rowSds(filtCounts),
     main = "Raw Counts sd vs mean",
     xlim = c(0,10000),
     ylim = c(0,5000))

# log2 transform

logCounts <- log2(filtCounts + 1)

statusCols <- str_replace_all(sampleinfo$Status, 
                              c(Infected = "red", Uninfected = "orange"))
boxplot(logCounts,
        xlab="",
        ylab="Log2(Counts)",
        las = 2,
        col=statusCols,
        main = "Log2(Counts)")

plot(rowMeans(logCounts), rowSds(logCounts),
     main = "Log2 Counts sd vs mean")

# vst

vstCounts <- vst(filtCounts)

boxplot(vstCounts,
        xlab="",
        ylab="vst(Counts)",
        las = 2,
        col=statusCols,
        main = "vst(Counts)")

plot(rowMeans(vstCounts), rowSds(vstCounts),
     main = "vst Counts sd vs mean")

# Exercise 2

rlogCounts <- rlog(filtCounts)

boxplot(rlogCounts,
        xlab="",
        ylab="rlog(Counts)",
        las = 2,
        col=statusCols,
        main = "rlog(Counts)")

# Principle Component Analysis

library(ggfortify)

pcDat <- prcomp(t(rlogCounts))
autoplot(pcDat)

autoplot(pcDat,
         data = sampleinfo,
         colour="Status",
         shape="TimePoint",
         size = 5)

# Exercise 3

autoplot(pcDat,
         data = sampleinfo,
         colour="Status",
         shape="TimePoint",
         size = 5,
         x = 2,
         y = 3)

# adding labels

library(ggrepel)

autoplot(pcDat,
         data = sampleinfo,  
         colour="Status", 
         shape="TimePoint",
         size=5) +
  geom_text_repel(aes(x=PC1, y=PC2, label=SampleName), box.padding = 0.8)

sampleinfo <- mutate(sampleinfo, Status=case_when(
  SampleName=="SRR7657882" ~ "Uninfected",
  SampleName=="SRR7657873" ~ "Infected", 
  TRUE ~ Status))


