rm(list = ls())

##change path that packages are loaded to
.libPaths("E:/R-packages") #make sure this is set up locally
.libPaths()#prints out the paths, in position 1 is where the packages will be loaded

options(lib.loc = "E:/R-packages")

## I don't need all these packages anymore, but too lazy to delete them so whatever.

library(ggplot2)
library(tibble)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(ggplot2)
library(tidybulk)
library(ComplexHeatmap)
library(tidyHeatmap)
library(ggrepel)
library(plotly)
library(RColorBrewer)
library(pheatmap)
library(apeglm)
library(ashr)
library(annotables)
library(tidyverse)
library(reshape2)
library(ggsignif)
library(statix)
library(ggpubr)
library(broom)
library(multcomp)

##upload logcpm
logcpm_counts <- read.csv('E:/paper-files/mlcm_total_logcpmc.csv', sep=',', header = TRUE)
colnames(logcpm_counts)[1] = "ensgene" #change column 1 name because imports wonky


###annotate logcpm+c file
## annotate significant genes
logcpm_counts <- data.frame(logcpm_counts) 
logcpm_counts <- left_join(x = logcpm_counts,
                          y = grcm38 [, c("ensgene", "symbol", "entrez", "biotype", "description")], 
                          by = "ensgene")

write.csv(logcpm_counts, "E:/paper-files/mlcm_total_logcpmc_annotated.csv", row.names = TRUE)

logcpm_counts <- filter(logcpm_counts, biotype == "protein_coding")

gc_markers <- c("Cyp19a1", "Fshr", "Tgfbr1", "Bmpr1b", "Foxl2", "Ahr", "Bmpr2", "Cyp11a1", "Star", "Inhba", 
                "Amh", "Kitl", "Hif1a")

logcpm_gcs <- filter(logcpm_counts, symbol %in% gc_markers)

epi_markers <- c("Cdh1", "Epcam", "Muc16", "Krt7", "Krt8", "Krt18", "Krt19", "Muc1", "Hsd17b10", 
                 "Cdh2", "Lgr5", "Ly6a", "Aldh1a1", "Wnt5a", "Cd44", "Aldh1a2", "Msln", "Egfr", "Pax8")

logcpm_epi <- filter(logcpm_counts, symbol %in% epi_markers)

condition <- c("Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma",
               "KO_stroma", "KO_stroma", "KO_stroma", "KO_stroma", "KO_stroma", "KO_stroma", "KO_stroma",
               "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords",
               "Mature GCs", "Mature GCs", "Mature GCs", "Mature GCs", "Mature GCs",  "Mature GCs",  "Mature GCs", 
               "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs", 
               "WT_stroma", "WT_stroma", "WT_stroma", "WT_stroma", "WT_stroma", "WT_stroma", "WT_stroma")

selected_columns <- logcpm_gcs[, c("Sxcd_1", "Sxcd_2", "Sxcd_3", "Sxcd_4", "Sxcd_5", "Sxcd_6", "Sxcd_7",
                                      "GC_1", "GC_2", "GC_3", "GC_4", "GC_5", "GC_6", "GC_7",
                                      "Prim_1", "Prim_2", "Prim_3", "Prim_4", "Prim_5", "Prim_6", "Prim_7")]

logcpm_gcs <- logcpm_gcs[, c(16:36, 44)]
logcpm_gcs$symbol <- str_to_title(logcpm_gcs$symbol)
logcpm_gcs$marker <- "GCs"
logcpm_epi <- logcpm_epi[, c(16:36, 44)]
logcpm_epi$symbol <- str_to_title(logcpm_epi$symbol)
logcpm_epi$marker <- "Epithelial"

logcpm_markers <- dplyr::bind_rows(logcpm_gcs, logcpm_epi)

###################################HEATMAP######################################
logcpm_markers <- unique(logcpm_markers)
rownames(logcpm_markers) <- logcpm_markers$symbol
logcpm_markers <- logcpm_markers[,1:21]
logcpm_markers_mat <- as.matrix(logcpm_markers)

sample_groups <- c("Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords",
                   "Mature GCs", "Mature GCs", "Mature GCs", "Mature GCs", "Mature GCs", "Mature GCs", "Mature GCs",
                   "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs")

sample_groups <- factor(sample_groups, levels = c("Sex cords", "Primary GCs", "Mature GCs"))

anno_colours <- c("Sex cords" = "#004D40", "Primary GCs" = "#5D286B", "Mature GCs" = "#1E88E5")

cat_colours <- c("GC Marker" = "grey50", "Epithelial Marker" = "springgreen")

categories <- c("GC Marker", "GC Marker", "GC Marker", "GC Marker", "GC Marker", 
                "GC Marker", "GC Marker", "GC Marker", "GC Marker", "GC Marker", 
                "GC Marker", "GC Marker", "GC Marker", 
                "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", 
                "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", 
                "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", 
                "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", "Epithelial Marker")

top_ha <-  HeatmapAnnotation(
  sample_groups = sample_groups, 
  col = list(sample_groups = anno_colours),
  annotation_height = unit(3, "mm"), 
  annotation_width = unit(0.5, 'cm'), 
  gap = unit(0.5, 'mm'), 
  border = TRUE, 
  show_annotation_name = FALSE, 
  annotation_legend_param = list(
    sample_groups = list(
      nrow = 3, 
      title = "Condition", 
      title_position = 'topleft',
      legend_direction = 'horizontal',
      title_gp = gpar(fontsize =12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'plain'))))

top_row <- ha_row <- HeatmapAnnotation(which = "row",
                                       categories = categories, 
                                       col = list(categories = cat_colours),
                                       annotation_height = 0.3, 
                                       annotation_width = unit(1, 'cm'), 
                                       gap = unit(1, 'mm'), 
                                       border = TRUE, 
                                       show_legend = TRUE, 
                                       show_annotation_name = FALSE, 
                                       annotation_legend_param = list(
                                         categories = list(
                                           nrow = 2,
                                           title = 'Type',
                                           title_position = 'topleft',
                                           legend_direction = 'vertical',
                                           title_gp = gpar(fontsize = 12, fontface = 'bold'),
                                           labels_gp = gpar(fontsize = 12, fontface = 'plain'))))

heatmap_markers <- ComplexHeatmap::Heatmap(logcpm_markers_mat, 
                                           cluster_rows = TRUE, 
                                           cluster_columns = TRUE, 
                                           name = "log2CPM + c", 
                                           row_km = 3, 
                                           rect_gp = gpar(col = "grey10", lwd = 0.25), 
                                           clustering_distance_rows = "euclidean", 
                                           row_names_gp = gpar(fontsize = 12, fontface = "italic"),
                                           show_row_names = TRUE, 
                                           show_column_names = FALSE, 
                                           heatmap_legend_param = list(title_gp = gpar(fontsize = 12,  fontface = "bold"), 
                                                                       labels_gp = gpar(fontsize = 12,  fontface = "plain")), 
                                           top_annotation = top_ha, 
                                           right_annotation = top_row)

#draw the heatmap
draw(heatmap_markers,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE)

                                           
