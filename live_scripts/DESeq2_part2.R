library(DESeq2)
library(tidyverse)

# recap of part 1
txi <- readRDS("RObjects/txi.rds")
sampleinfo <- read_tsv("data/samplesheet_corrected.tsv", col_types = "cccc") %>%
  mutate(Status = fct_relevel(Status, "Uninfected"))
simple.model <- as.formula(~ Status)
ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = simple.model)
keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep,]

ddsObj <- DESeq(ddsObj.filt)

results.simple <- results(ddsObj, alpha = 0.05)

results.simple


# Exercise 2

additive.model <- as.formula(~ TimePoint + Status)

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = additive.model)
keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep,]

ddsObj <- DESeq(ddsObj.filt)

results.additive <- results(ddsObj, alpha = 0.05)
results.additive

# why do we get this contrast?

model.matrix(additive.model, data = sampleinfo)

resultsNames(ddsObj)

results.InfectedvUninfected <- results.additive
rm(results.additive)

# what are our top sig diff genes?

topGenesIvU <- as.data.frame(results.InfectedvUninfected) %>%
  rownames_to_column("GeneID") %>%
  top_n(100, wt = -padj) %>%
  arrange(padj)
  
# Exercise 3

resultsNames(ddsObj)
results.d33vd11 <- results(ddsObj, name = "TimePoint_d33_vs_d11", alpha = 0.05)
results.d33vd11
  
sum(results.d33vd11$padj < 0.05, na.rm = TRUE)

# lets look at our PCA? Do we want to try different model?

vstcounts <- vst(ddsObj.raw, blind = TRUE)
plotPCA(vstcounts, intgroup = c("Status", "TimePoint"))

# How do we test between models?

ddsObj.LRT <- DESeq(ddsObj, test = "LRT", reduced = simple.model)

results.Additive_v_Simple <- results(ddsObj.LRT, alpha = 0.05)
results.Additive_v_Simple

sum(results.Additive_v_Simple$padj < 0.05, na.rm = TRUE)

# Exercise 4

interaction.model <- as.formula(~TimePoint * Status)

ddsObj.raw <- DESeqDataSetFromTximport(txi = txi,
                                       colData = sampleinfo,
                                       design = interaction.model)
keep <- rowSums(counts(ddsObj.raw)) > 5
ddsObj.filt <- ddsObj.raw[keep,]

ddsObj.interaction <- DESeq(ddsObj.filt)

# There are lots more possible contrasts, how do I access them?

resultsNames(ddsObj.interaction)

results.interaction.11 <- results(ddsObj.interaction, 
                                  name = "Status_Infected_vs_Uninfected",
                                  alpha = 0.05) 
results.interaction.11

results.interaction.33 <- results(ddsObj.interaction,
                                  contrast = list(c("Status_Infected_vs_Uninfected",
                                                    "TimePointd33.StatusInfected")),
                                  alpha = 0.05)  
results.interaction.33

sum(results.interaction.11$padj < 0.05, na.rm = TRUE)
sum(results.interaction.33$padj < 0.05, na.rm = TRUE)

# testing interaction vs. additive

ddsObj.interaction.LRT <- DESeq(ddsObj.interaction, test = "LRT", reduced = additive.model)

results.interaction_v_additive <- results(ddsObj.interaction.LRT, alpha = 0.05)


sum(results.interaction_v_additive$padj < 0.05, na.rm = TRUE)

# Exercise 5

results.d33_v_d11_uninfected <- results(ddsObj.interaction, 
                                        name = "TimePoint_d33_vs_d11",
                                        alpha = 0.05)
table(results.d33_v_d11_uninfected$padj < 0.05)

results.d33_v_d11_infected <- results(ddsObj.interaction,
                                      contrast = list(c("TimePoint_d33_vs_d11", "TimePointd33.StatusInfected")),
                                      alpha = 0.05)

table(results.d33_v_d11_infected$padj < 0.05)






