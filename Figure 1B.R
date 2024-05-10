###################### Figure 1B

## Data preparation
pval_threshold <- 0.05

Gapmer_DEGS <- rownames(results_GapmervsUT.df[which(results_GapmervsUT.df$padj < pval_threshold), ])
Gapmer_DEGS

PMO_DEGS <- rownames(results_PMOvsUT.df[which(results_PMOvsUT.df$padj < pval_threshold), ])
PMO_DEGS

TMO_DEGS <- rownames(results_TMOvsUT.df[which(results_TMOvsUT.df$padj < pval_threshold), ])
TMO_DEGS

MOE_DEGS <- rownames(results_MOEvsUT.df[which(results_MOEvsUT.df$padj < pval_threshold), ])
MOE_DEGS



GapmerControl_DEGS <- rownames(results_GapmerControlvsUT.df[which(results_GapmerControlvsUT.df$padj < pval_threshold), ])
GapmerControl_DEGS

PMO_GTC_DEGS <- rownames(results_PMO_GTCvsUT.df[which(results_PMO_GTCvsUT.df$padj < pval_threshold), ])
PMO_GTC_DEGS

TMO_GTC_DEGS <- rownames(results_TMO_GTCvsUT.df[which(results_TMO_GTCvsUT.df$padj < pval_threshold), ])
TMO_GTC_DEGS

MOE_GTC_DEGS <- rownames(results_MOE_GTCvsUT.df[which(results_MOE_GTCvsUT.df$padj < pval_threshold), ])
MOE_GTC_DEGS

################# Figure 1B

library(VennDiagram)

custom_fill <- c("salmon1", "palegreen", "lightskyblue", "plum1")


venn.plot1 <- venn.diagram(
  x = list(Gapmer = Gapmer_DEGS, PMO = PMO_DEGS, TMO = TMO_DEGS, MOE = MOE_DEGS),
  category.names = c("Gapmer", "PMO", "TMO", "MOE"),
  filename = NULL,
  fill = custom_fill
  
)

grid.draw(venn.plot1)

########################## Supplementary Figure Upset plot
library(UpSetR)
library(ComplexHeatmap)

# Define your list of gene sets
x <- list(
  Gapmer_ITGA4 = Gapmer_DEGS,
  Gapmer_Control = GapmerControl_DEGS,
  PMO_ITGA4 = PMO_DEGS,
  PMO_GTC = PMO_GTC_DEGS,
  TMO_ITGA4 = TMO_DEGS,
  TMO_GTC = TMO_GTC_DEGS,
  MOE_ITGA4 = MOE_DEGS,
  MOE_GTC = MOE_GTC_DEGS
)

# Create the combination matrix
m <- make_comb_mat(x)


set_order <- c("Gapmer_ITGA4", "Gapmer_Control", "PMO_ITGA4", "PMO_GTC", 
               "TMO_ITGA4", "TMO_GTC", "MOE_ITGA4", "MOE_GTC")

ss = set_size(m)
cs = comb_size(m)

ht = UpSet(m, 
           set_order = set_order,
           comb_order = order(comb_degree(m), -cs),
           top_annotation = HeatmapAnnotation(
             "Intersection Size" = anno_barplot(cs, 
                                                ylim = c(0, max(cs)*1.1),
                                                border = FALSE, 
                                                gp = gpar(fill = "black"), 
                                                height = unit(4, "cm")
             ), 
             annotation_name_side = "left", 
             annotation_name_rot = 90),
           left_annotation = rowAnnotation(
             "Set Size" = anno_barplot(-ss, 
                                       baseline = 0,
                                       axis_param = list(
                                         at = c(0, -1000, -2000, -3000, -4000),
                                         labels = c(0, 1000, 2000, 3000, 4000),
                                         labels_rot = 0),
                                       border = FALSE, 
                                       gp = gpar(fill = "black"), 
                                       width = unit(4, "cm")
             ),
             set_name = anno_text(set_name(m), 
                                  location = 0.5, 
                                  just = "center",
                                  width = max_text_width(set_name(m)) + unit(4, "mm"))
           ), 
           right_annotation = NULL,
           show_row_names = FALSE)
ht = draw(ht)
od = column_order(ht)
decorate_annotation("Intersection Size", {
  grid.text(cs[od], x = seq_along(cs), y = unit(cs[od], "native") + unit(2, "pt"), 
            default.units = "native", just = c("center", "bottom"), 
            gp = gpar(fontsize = 8, col = "black"), rot = 0)
})

