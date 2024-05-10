library(openxlsx)

# File paths
RawCounts <- "Data/RawCounts.xlsx"
SampleInfo <- "Data/SampleInfo.xlsx"

# Import counts data
RawCounts_sheets <- openxlsx::getSheetNames(RawCounts)

RawCounts_df <- lapply(RawCounts_sheets, function(sheet) {
  counts <- openxlsx::read.xlsx(RawCounts, sheet)
  row.names(counts) <- counts$GeneID
  counts$GeneID <- NULL
  counts
})
names(RawCounts_df) <- RawCounts_sheets

# Import sample info 
SampleInfo_sheets <- openxlsx::getSheetNames(SampleInfo)

SampleInfo_df <- lapply(SampleInfo_sheets, function(sheet) {
  info <- openxlsx::read.xlsx(SampleInfo, sheet)
  row.names(info) <- info$Sample
  info$Sample <- NULL
  info
})
names(SampleInfo_df) <- SampleInfo_sheets


### Differential gene expression 
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)


dds_GapmervsUT <- DESeqDataSetFromMatrix(countData = RawCounts_df$GapmervsUT,
                                         colData = SampleInfo_df$GapmervsUT_info,
                                         design =  ~ Condition)
dds_GapmerControlvsUT <- DESeqDataSetFromMatrix(countData = RawCounts_df$GapmerControlvsUT,
                                                 colData = SampleInfo_df$GapmerControlvsUT_info,
                                                 design =  ~ Condition)

dds_PMOvsUT <- DESeqDataSetFromMatrix(countData = RawCounts_df$PMOvsUT,
                                     colData = SampleInfo_df$PMOvsUT_info,
                                     design =  ~ Condition)
dds_PMO_GTCvsUT <- DESeqDataSetFromMatrix(countData = RawCounts_df$PMO_GTCvsUT,
                                          colData = SampleInfo_df$PMO_GTCvsUT_info,
                                          design =  ~ Condition)

dds_TMOvsUT <- DESeqDataSetFromMatrix(countData = RawCounts_df$TMOvsUT,
                                      colData = SampleInfo_df$TMOvsUT_info,
                                      design =  ~ Condition)
dds_TMO_GTCvsUT <- DESeqDataSetFromMatrix(countData = RawCounts_df$TMO_GTCvsUT,
                                          colData = SampleInfo_df$TMO_GTCvsUT_info,
                                          design =  ~ Condition)

dds_MOEvsUT <- DESeqDataSetFromMatrix(countData = RawCounts_df$MOEvsUT,
                                      colData = SampleInfo_df$MOEvsUT_info,
                                      design =  ~ Condition)
dds_MOE_GTCvsUT <- DESeqDataSetFromMatrix(countData = RawCounts_df$MOE_GTCvsUT,
                                          colData = SampleInfo_df$MOE_GTCvsUT_info,
                                          design =  ~ Condition)

#################################################################

keep_GapmervsUT <- rowSums(counts(dds_GapmervsUT) >= 10) >=3
dds_GapmervsUT <- dds_GapmervsUT[keep_GapmervsUT,]

keep_GapmerControlvsUT <- rowSums(counts(dds_GapmerControlvsUT) >= 10) >=3
dds_GapmerControlvsUT <- dds_GapmerControlvsUT[keep_GapmerControlvsUT,]

keep_PMOvsUT <- rowSums(counts(dds_PMOvsUT) >= 10) >=3
dds_PMOvsUT <- dds_PMOvsUT[keep_PMOvsUT,]

keep_PMO_GTCvsUT <- rowSums(counts(dds_PMO_GTCvsUT) >= 10) >=3
dds_PMO_GTCvsUT <- dds_PMO_GTCvsUT[keep_PMO_GTCvsUT,]

keep_TMOvsUT <- rowSums(counts(dds_TMOvsUT) >= 10) >=3
dds_TMOvsUT <- dds_TMOvsUT[keep_TMOvsUT,]

keep_TMO_GTCvsUT <- rowSums(counts(dds_TMO_GTCvsUT) >= 10) >=3
dds_TMO_GTCvsUT <- dds_TMO_GTCvsUT[keep_TMO_GTCvsUT,]

keep_MOEvsUT <- rowSums(counts(dds_MOEvsUT) >= 10) >=3
dds_MOEvsUT <- dds_MOEvsUT[keep_MOEvsUT,]

keep_MOE_GTCvsUT <- rowSums(counts(dds_MOE_GTCvsUT) >= 10) >=3
dds_MOE_GTCvsUT <- dds_MOE_GTCvsUT[keep_MOE_GTCvsUT,]

####################################################

deseq_GapmervsUT <- DESeq(dds_GapmervsUT)
deseq_GapmerControlvsUT <- DESeq(dds_GapmerControlvsUT)

deseq_PMOvsUT <- DESeq(dds_PMOvsUT)
deseq_PMO_GTCvsUT <- DESeq(dds_PMO_GTCvsUT)

deseq_TMOvsUT <- DESeq(dds_TMOvsUT)
deseq_TMO_GTCvsUT <- DESeq(dds_TMO_GTCvsUT)

deseq_MOEvsUT <- DESeq(dds_MOEvsUT)
deseq_MOE_GTCvsUT <- DESeq(dds_MOE_GTCvsUT)
#######################################################
results_GapmervsUT <- results(deseq_GapmervsUT, pAdjustMethod="BH", contrast = c("Condition", "Gapmer", "Untreated"))
results_GapmerControlvsUT <- results(deseq_GapmerControlvsUT, pAdjustMethod="BH", contrast = c("Condition", "GapmerControl", "Untreated"))

results_PMOvsUT <- results(deseq_PMOvsUT, pAdjustMethod="BH", contrast = c("Condition", "PMO", "Untreated"))
results_PMO_GTCvsUT <- results(deseq_PMO_GTCvsUT, pAdjustMethod="BH", contrast = c("Condition", "PMO_GTC", "Untreated"))

results_TMOvsUT <- results(deseq_TMOvsUT, pAdjustMethod="BH", contrast = c("Condition", "TMO", "Untreated"))
results_TMO_GTCvsUT <- results(deseq_TMO_GTCvsUT, pAdjustMethod="BH", contrast = c("Condition", "TMO_GTC", "Untreated"))

results_MOEvsUT <- results(deseq_MOEvsUT, pAdjustMethod="BH", contrast = c("Condition", "MOE", "Untreated"))
results_MOE_GTCvsUT <- results(deseq_MOE_GTCvsUT, pAdjustMethod="BH", contrast = c("Condition", "MOE_GTC", "Untreated"))

###########################################################
#Sort by p-value

results_GapmervsUT <- results_GapmervsUT[order(results_GapmervsUT$padj),]
results_GapmerControlvsUT <- results_GapmerControlvsUT[order(results_GapmerControlvsUT$padj),]

results_PMOvsUT <- results_PMOvsUT[order(results_PMOvsUT$padj),]
results_PMO_GTCvsUT <- results_PMO_GTCvsUT[order(results_PMO_GTCvsUT$padj),]

results_TMOvsUT <- results_TMOvsUT[order(results_TMOvsUT$padj),]
results_TMO_GTCvsUT <- results_TMO_GTCvsUT[order(results_TMO_GTCvsUT$padj),]

results_MOEvsUT <- results_MOEvsUT[order(results_MOEvsUT$padj),]
results_MOE_GTCvsUT <- results_MOE_GTCvsUT[order(results_MOE_GTCvsUT$padj),]
##############################################      

results_GapmervsUT.df <- as.data.frame(results_GapmervsUT)
results_GapmerControlvsUT.df <- as.data.frame(results_GapmerControlvsUT)

results_PMOvsUT.df <- as.data.frame(results_PMOvsUT)
results_PMO_GTCvsUT.df <- as.data.frame(results_PMO_GTCvsUT)

results_TMOvsUT.df <- as.data.frame(results_TMOvsUT)
results_TMO_GTCvsUT.df <- as.data.frame(results_TMO_GTCvsUT)

results_MOEvsUT.df <- as.data.frame(results_MOEvsUT)
results_MOE_GTCvsUT.df <- as.data.frame(results_MOE_GTCvsUT)
#################################################
write.csv(results_GapmervsUT.df, file = "Data/DESeq2_GapmervsUT.csv")
write.csv(results_GapmerControlvsUT.df, file = "Data/DESeq2_GapmerControlvsUT.csv")

write.csv(results_PMOvsUT.df, file = "Data/DESeq2_PMOvsUT.csv")
write.csv(results_PMO_GTCvsUT.df, file = "Data/DESeq2_PMO_GTCvsUT.csv")

write.csv(results_TMOvsUT.df, file = "Data/DESeq2_TMOvsUT.csv")
write.csv(results_TMO_GTCvsUT.df, file = "Data/DESeq2_TMO_GTCvsUT.csv")

write.csv(results_MOEvsUT.df, file = "Data/DESeq2_MOEvsUT.csv")
write.csv(results_MOE_GTCvsUT.df, file = "Data/DESeq2_MOE_GTCvsUT.csv")

######################################### Figure 1A

GapmervsUT_volc <- EnhancedVolcano(
  results_GapmervsUT.df,
  lab = rownames(results_GapmervsUT.df),
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-10, 10),
  ylim = c(0, 20),
  pCutoff = 0.05,
  FCcutoff = 0.2,
  col = c("black", "black", "black", "red"),
  colAlpha = 1,
  labSize = 3,
  legendPosition = "",
  legendLabSize = 10,
  caption = NULL,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  title = "Gapmer",
  subtitle = "",
  border = 'full',
  borderWidth = 0.8,
  borderColour = 'black')


GapmervsUT_volc <- GapmervsUT_volc + 
  xlab("log2FoldChange") +
  ylab("-log10(padj)") +
  theme(plot.title = element_text(margin = margin(b = -18)))

print(GapmervsUT_volc)

###

PMOvsUT_volc <- EnhancedVolcano(
  results_PMOvsUT.df,
  lab = rownames(results_PMOvsUT.df),
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-10, 10),
  ylim = c(0, 20),
  pCutoff = 0.05,
  FCcutoff = 0.2,
  col = c("black", "black", "black", "red"),
  colAlpha = 1,
  labSize = 3,
  legendPosition = "",
  legendLabSize = 10,
  caption = NULL,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  title = "PMO",
  subtitle = "",
  border = 'full',
  borderWidth = 0.8,
  borderColour = 'black')


PMOvsUT_volc <- PMOvsUT_volc + 
  xlab("log2FoldChange") +
  ylab("-log10(padj)") +
  theme(plot.title = element_text(margin = margin(b = -18)))

print(PMOvsUT_volc)
##
TMOvsUT_volc <- EnhancedVolcano(
  results_TMOvsUT.df,
  lab = rownames(results_TMOvsUT.df),
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-10, 10),
  ylim = c(0, 20),
  pCutoff = 0.05,
  FCcutoff = 0.2,
  col = c("black", "black", "black", "red"),
  colAlpha = 1,
  labSize = 3,
  legendPosition = "",
  legendLabSize = 10,
  caption = NULL,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  title = "TMO",
  subtitle = "",
  border = 'full',
  borderWidth = 0.8,
  borderColour = 'black')


TMOvsUT_volc <- TMOvsUT_volc + 
  xlab("log2FoldChange") +
  ylab("-log10(padj)") +
  theme(plot.title = element_text(margin = margin(b = -18)))

print(TMOvsUT_volc)
###
MOEvsUT_volc <- EnhancedVolcano(
  results_MOEvsUT.df,
  lab = rownames(results_MOEvsUT.df),
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-10, 10),
  ylim = c(0, 20),
  pCutoff = 0.05,
  FCcutoff = 0.2,
  col = c("black", "black", "black", "red"),
  colAlpha = 1,
  labSize = 3,
  legendPosition = "",
  legendLabSize = 10,
  caption = NULL,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  title = "MOE",
  subtitle = "",
  border = 'full',
  borderWidth = 0.8,
  borderColour = 'black')


MOEvsUT_volc <- MOEvsUT_volc + 
  xlab("log2FoldChange") +
  ylab("-log10(padj)") +
  theme(plot.title = element_text(margin = margin(b = -18)))

print(MOEvsUT_volc)

######################### Supplementary GTC
GapmerControlvsUT_volc <- EnhancedVolcano(
  results_GapmerControlvsUT.df,
  lab = rownames(results_GapmerControlvsUT.df),
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-10, 10),
  ylim = c(0, 20),
  pCutoff = 0.05,
  FCcutoff = 0.2,
  col = c("black", "black", "black", "red"),
  colAlpha = 1,
  labSize = 3,
  legendPosition = "",
  legendLabSize = 10,
  caption = NULL,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  title = "Gapmer Control",
  subtitle = "",
  border = 'full',
  borderWidth = 0.8,
  borderColour = 'black')

GapmerControlvsUT_volc <- GapmerControlvsUT_volc + 
  xlab("log2FoldChange") +
  ylab("-log10(padj)") +
  theme(plot.title = element_text(margin = margin(b = -18))) 

print(GapmerControlvsUT_volc)
###
PMO_GTCvsUT_volc <- EnhancedVolcano(
  results_PMO_GTCvsUT.df,
  lab = rownames(results_PMO_GTCvsUT.df),
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-10, 10),
  ylim = c(0, 20),
  pCutoff = 0.05,
  FCcutoff = 0.2,
  col = c("black", "black", "black", "red"),
  colAlpha = 1,
  labSize = 3,
  legendPosition = "",
  legendLabSize = 10,
  caption = NULL,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  title = "PMO GTC",
  subtitle = "",
  border = 'full',
  borderWidth = 0.8,
  borderColour = 'black')


PMO_GTCvsUT_volc <- PMO_GTCvsUT_volc + 
  xlab("log2FoldChange") +
  ylab("-log10(padj)") +
  theme(plot.title = element_text(margin = margin(b = -18)))

print(PMO_GTCvsUT_volc)

####
TMO_GTCvsUT_volc <- EnhancedVolcano(
  results_TMO_GTCvsUT.df,
  lab = rownames(results_TMO_GTCvsUT.df),
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-10, 10),
  ylim = c(0, 20),
  pCutoff = 0.05,
  FCcutoff = 0.2,
  col = c("black", "black", "black", "red"),
  colAlpha = 1,
  labSize = 3,
  legendPosition = "",
  legendLabSize = 10,
  caption = NULL,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  title = "TMO GTC",
  subtitle = "",
  border = 'full',
  borderWidth = 0.8,
  borderColour = 'black')


TMO_GTCvsUT_volc <- TMO_GTCvsUT_volc + 
  xlab("log2FoldChange") +
  ylab("-log10(padj)") +
  theme(plot.title = element_text(margin = margin(b = -18)))

print(TMO_GTCvsUT_volc)
##
MOE_GTCvsUT_volc <- EnhancedVolcano(
  results_MOE_GTCvsUT.df,
  lab = rownames(results_MOE_GTCvsUT.df),
  x = 'log2FoldChange',
  y = 'padj',
  xlim = c(-10, 10),
  ylim = c(0, 20),
  pCutoff = 0.05,
  FCcutoff = 0.2,
  col = c("black", "black", "black", "red"),
  colAlpha = 1,
  labSize = 3,
  legendPosition = "",
  legendLabSize = 10,
  caption = NULL,
  gridlines.major = TRUE,
  gridlines.minor = TRUE,
  title = "MOE GTC",
  subtitle = "",
  border = 'full',
  borderWidth = 0.8,
  borderColour = 'black')


MOE_GTCvsUT_volc <- MOE_GTCvsUT_volc + 
  xlab("log2FoldChange") +
  ylab("-log10(padj)") +
  theme(plot.title = element_text(margin = margin(b = -18)))

print(MOE_GTCvsUT_volc)

#######################################################

library(ggpubr)
ggarrange(GapmervsUT_volc, PMOvsUT_volc, TMOvsUT_volc, MOEvsUT_volc, ncol = 4, nrow = 1)