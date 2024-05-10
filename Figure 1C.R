########################## Figure 1C
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)


sig_GapmervsUT <- na.omit(results_GapmervsUT.df)
sig_GapmervsUT <- sig_GapmervsUT[sig_GapmervsUT$padj < 0.05 & abs(sig_GapmervsUT$log2FoldChange) >= 0.2,]

sig_GapmerControlvsUT <- na.omit(results_GapmerControlvsUT.df)
sig_GapmerControlvsUT <- sig_GapmerControlvsUT[sig_GapmerControlvsUT$padj < 0.05 & abs(sig_GapmerControlvsUT$log2FoldChange) >= 0.2,]

sig_PMO_GTCvsUT <- na.omit(results_PMO_GTCvsUT.df)
sig_PMO_GTCvsUT <- sig_PMO_GTCvsUT[sig_PMO_GTCvsUT$padj < 0.05 & abs(sig_PMO_GTCvsUT$log2FoldChange) >= 0.2,]

sig_MOEvsUT <- na.omit(results_MOEvsUT.df)
sig_MOEvsUT <- sig_MOEvsUT[sig_MOEvsUT$padj < 0.05 & abs(sig_MOEvsUT$log2FoldChange) >= 0.2,]

sig_MOE_GTCvsUT <- na.omit(results_MOE_GTCvsUT.df)
sig_MOE_GTCvsUT <- sig_MOE_GTCvsUT[sig_MOE_GTCvsUT$padj < 0.05 & abs(sig_MOE_GTCvsUT$log2FoldChange) >= 0.2,]


####################################################


#Getting log2 fold change
Gapmer_gene_list <- sig_GapmervsUT$log2FoldChange
GapmerControl_gene_list <- sig_GapmerControlvsUT$log2FoldChange
PMO_GTC_gene_list <- sig_PMO_GTCvsUT$log2FoldChange
MOE_gene_list <- sig_MOEvsUT$log2FoldChange
MOE_GTC_gene_list <- sig_MOE_GTCvsUT$log2FoldChange

# name the vector
names(Gapmer_gene_list) <- rownames(sig_GapmervsUT)
names(GapmerControl_gene_list) <- rownames(sig_GapmerControlvsUT)
names(PMO_GTC_gene_list) <- rownames(sig_PMO_GTCvsUT)
names(MOE_gene_list) <- rownames(sig_MOEvsUT)
names(MOE_GTC_gene_list) <- rownames(sig_MOE_GTCvsUT)



# omit any NA values 
Gapmer_gene_list<-na.omit(Gapmer_gene_list)
GapmerControl_gene_list<-na.omit(GapmerControl_gene_list)
PMO_GTC_gene_list<-na.omit(PMO_GTC_gene_list)
MOE_gene_list<-na.omit(MOE_gene_list)
MOE_GTC_gene_list<-na.omit(MOE_GTC_gene_list)

# sort the list in decreasing order
Gapmer_gene_list = sort(Gapmer_gene_list, decreasing = TRUE)
GapmerControl_gene_list = sort(GapmerControl_gene_list, decreasing = TRUE)
PMO_GTC_gene_list = sort(PMO_GTC_gene_list, decreasing = TRUE)
MOE_gene_list = sort(MOE_gene_list, decreasing = TRUE)
MOE_GTC_gene_list = sort(MOE_GTC_gene_list, decreasing = TRUE)

##############################################################
set.seed(123)

gse_Gapmer <- gseGO(geneList=Gapmer_gene_list, 
              ont ="ALL", 
              keyType = "SYMBOL", 
              minGSSize = 3, 
              maxGSSize = 800, 
              pvalueCutoff = 0.05, 
              verbose = TRUE, 
              OrgDb = "org.Hs.eg.db")  
              

gse_GapmerControl <- gseGO(geneList=GapmerControl_gene_list, 
              ont ="ALL", 
              keyType = "SYMBOL", 
              minGSSize = 3, 
              maxGSSize = 800, 
              pvalueCutoff = 0.05, 
              verbose = TRUE, 
              OrgDb = "org.Hs.eg.db")

#NOTE : no term enriched under specific pvalueCutoff (0.05) for GapmerControl

gse_PMO_GTC <- gseGO(geneList=PMO_GTC_gene_list, 
                     ont ="ALL", 
                     keyType = "SYMBOL", 
                     minGSSize = 3, 
                     maxGSSize = 800, 
                     pvalueCutoff = 0.05, 
                     verbose = TRUE, 
                     OrgDb = "org.Hs.eg.db")

#NOTE : no term enriched under specific pvalueCutoff (0.05) for PMO_GTC

gse_MOE <- gseGO(geneList=MOE_gene_list, 
              ont ="ALL", 
              keyType = "SYMBOL", 
              minGSSize = 3, 
              maxGSSize = 800, 
              pvalueCutoff = 0.05, 
              verbose = TRUE, 
              OrgDb = "org.Hs.eg.db")

gse_MOE_GTC <- gseGO(geneList=MOE_GTC_gene_list, 
              ont ="ALL", 
              keyType = "SYMBOL", 
              minGSSize = 3, 
              maxGSSize = 800, 
              pvalueCutoff = 0.05, 
              verbose = TRUE, 
              OrgDb = "org.Hs.eg.db")


write.csv(gse_Gapmer@result, file = "GSE_Gapmer.csv")
write.csv(gse_GapmerControl@result, file = "GSE_GapmerControl.csv")
write.csv(gse_PMO_GTC@result, file = "GSE_PMO_GTC.csv")
write.csv(gse_MOE@result, file = "GSE_MOE.csv")
write.csv(gse_MOE_GTC@result, file = "GSE_MOE_GTC.csv")

############################################################## Figure 1C

require(DOSE)

dotplot(gse_Gapmer, showCategory = 15, font.size = 10, label_format = 40, title = "" , split=".sign") + 
  facet_grid(.~.sign)+ scale_x_continuous(name = "GeneRatio", limits = c(0, 1), breaks = seq(0.25, 1, by = 0.25)) +
  theme(plot.margin = margin(5, 5, 5, 5, "mm"))

dotplot(gse_MOE,showCategory = 15, font.size = 10, label_format = 40, title = "" , split=".sign") + 
  facet_grid(.~.sign)+ scale_x_continuous(name = "GeneRatio", limits = c(0, 1), breaks = seq(0.25, 1, by = 0.25)) +
  theme(plot.margin = margin(5, 5, 5, 5, "mm"))

##################################### Supplementary Figures

dotplot(gse_MOE_GTC, showCategory = 15, font.size = 10, label_format = 40, title = "" , split=".sign") + 
  facet_grid(.~.sign)+ scale_x_continuous(name = "GeneRatio", limits = c(0, 1), breaks = seq(0.25, 1, by = 0.25)) +
  theme(plot.margin = margin(5, 5, 5, 5, "mm"))

########################################################## Supplementary Figures

cnetplot(gse_Gapmer, categorySize="padjust", foldChange=Gapmer_gene_list, showCategory = 10)
cnetplot(gse_MOE, categorySize="padjust", foldChange=MOE_gene_list, showCategory = 10)
cnetplot(gse_MOE_GTC, categorySize="padjust", foldChange=MOE_GTC_gene_list, showCategory = 10)

