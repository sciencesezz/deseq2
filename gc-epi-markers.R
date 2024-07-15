rm(list = ls())

##change path that packages are loaded to
.libPaths("E:/R-packages") #make sure this is set up locally
.libPaths()#prints out the paths, in position 1 is where the packages will be loaded

options(lib.loc = "E:/R-packages")

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
logcpm_counts <- read.csv('E:/R-output/logcpm_mouse_mlcm.csv', sep=',', header = TRUE)
logcpm_counts <- filter(logcpm_counts, biotype == "protein_coding")

gc_markers <- c("cyp19a1", "fshr", "tgfbr1", "bmpr1b", "foxl2", "ahr", "bmpr2", "cyp11a1", "star", "inhba", 
                "amh", "kitl", "hif1a")

logcpm_gcs <- filter(logcpm_counts, symbol %in% gc_markers)

epi_markers <- c("cdh1", "epcam", "muc16", "krt7", "krt8", "krt18", "krt19", "muc1", "hsd17b10", 
                 "cdh2", "lgr5", "ly6a", "aldh1a1", "wnt5a", "cd44", "aldh1a2", "msln", "egfr", "pax8")

logcpm_epi <- filter(logcpm_counts, symbol %in% epi_markers)

condition <- c("adenoma", "adenoma", "adenoma", "adenoma", "adenoma", "adenoma", "adenoma",
               "KO_stroma", "KO_stroma", "KO_stroma", "KO_stroma", "KO_stroma", "KO_stroma", "KO_stroma",
               "sex cords", "sex cords", "sex cords", "sex cords", "sex cords", "sex cords", "sex cords",
               "mature GCs", "mature GCs", "mature GCs", "mature GCs", "mature GCs",  "mature GCs",  "mature GCs", 
               "primary GCs", "primary GCs", "primary GCs", "primary GCs", "primary GCs", "primary GCs", "primary GCs", 
               "WT_stroma", "WT_stroma", "WT_stroma", "WT_stroma", "WT_stroma", "WT_stroma", "WT_stroma")

selected_columns <- logcpm_gcs[, c("Sxcd_1", "Sxcd_2", "Sxcd_3", "Sxcd_4", "Sxcd_5", "Sxcd_6", "Sxcd_7",
                                      "GC_1", "GC_2", "GC_3", "GC_4", "GC_5", "GC_6", "GC_7",
                                      "Prim_1", "Prim_2", "Prim_3", "Prim_4", "Prim_5", "Prim_6", "Prim_7")]

logcpm_gcs <- logcpm_gcs[, c(16:36, 44)]
logcpm_gcs$symbol <- str_to_title(logcpm_gcs$symbol)
logcpm_gcs$marker <- "gcs"
logcpm_epi <- logcpm_epi[, c(16:36, 44)]
logcpm_epi$symbol <- str_to_title(logcpm_epi$symbol)
logcpm_epi$marker <- "epithelial"

logcpm_markers <- dplyr::bind_rows(logcpm_gcs, logcpm_epi)

#rename_cols <- c("Sxcds", "Sxcds", "Sxcds", "Sxcds", "Sxcds", "Sxcds", "Sxcds", 
                # "mature GCs", "mature GCs", "mature GCs", "mature GCs", "mature GCs", "mature GCs", "mature GCs", 
                # "primary GCs", "primary GCs", "primary GCs", "primary GCs", "primary GCs", "primary GCs", "primary GCs", 
                # "symbol", "marker")

#colnames(logcpm_markers) <- rename_cols

logcpm_gcs <- reshape2::melt(logcpm_gcs, id.vars = "symbol")
logcpm_epi <- reshape2::melt(logcpm_epi, id.vars = "symbol")

logcpm_gcs$value <- as.numeric(logcpm_gcs$value)
logcpm_epi$value <- as.numeric(logcpm_epi$value)

logcpm_gcs <- logcpm_gcs %>%
  mutate(
    variable = as.character(variable), 
    group = case_when(
    startsWith(variable, "Sxcd") ~ "Sex cords",
    startsWith(variable, "Prim") ~ "primary GCs",
    startsWith(variable, "GC") ~ "mature GCs",
    TRUE ~ "Other"
  ))

logcpm_epi <- logcpm_epi %>%
  mutate(
    variable = as.character(variable), 
    group = case_when(
      startsWith(variable, "Sxcd") ~ "Sex cords",
      startsWith(variable, "Prim") ~ "primary GCs",
      startsWith(variable, "GC") ~ "mature GCs",
      TRUE ~ "Other"
    ))

logcpm_gcs <- logcpm_gcs[1:273, ]
logcpm_epi <- logcpm_epi[1:420, ]
logcpm_gcs$group <- factor(logcpm_gcs$group, levels = c("Sex cords", "primary GCs", "mature GCs"))
logcpm_epi$group <- factor(logcpm_epi$group, levels = c("Sex cords", "primary GCs", "mature GCs"))

view(logcpm_gcs)
view(logcpm_epi)

gcs <- ggplot(logcpm_gcs, aes(x = group, y = value, colour = group)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.15, alpha = 0.5) + 
  facet_wrap(~ symbol, scales = "free_y") + 
  labs(x = "", y = "log2(CPM + c)") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank()) +
  theme(text = element_text(size = 16.0)) +
  theme(legend.position = "bottom")

print(gcs)
        
ggsave("E:/R-output/gc_markers.png", plot = gcs, height = 10, width = 8, dpi = 1200)

epi <- ggplot(logcpm_epi, aes(x = group, y = value, colour = group)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.15, alpha = 0.5) + 
  facet_wrap(~ symbol, scales = "free_y") + 
  labs(x = "", y = "log2(CPM + c)") + 
  theme_bw() + 
  theme(panel.grid = element_blank()) + 
  theme(axis.line.x = element_blank(),
        axis.text.x = element_blank()) +
  theme(text = element_text(size = 16.0)) +
  theme(legend.position = "bottom")

ggsave("E:/R-output/epi_markers.png", plot = epi, height = 10, width = 8, dpi = 1200)

###################################HEATMAP######################################
logcpm_markers <- unique(logcpm_markers)
rownames(logcpm_markers) <- logcpm_markers$symbol
logcpm_markers <- logcpm_markers[,1:21]
logcpm_markers_mat <- as.matrix(logcpm_markers)

sample_groups <- c("Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", 
                   "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs", 
                   "Mature GCs", "Mature GCs", "Mature GCs", "Mature GCs", "Mature GCs", "Mature GCs", "Mature GCs")

sample_groups <- factor(sample_groups, levels = c("Sex cords", "Primary GCs", "Mature GCs"))

anno_colours <- c("Sex cords" = "plum3", "Primary GCs" = "lightgreen", "Mature GCs" = "cadetblue")

sample_groups <- factor(sample_groups, levels = names(anno_colours))

categories <- c("GC Marker", "GC Marker", "GC Marker", "GC Marker", "GC Marker", 
                "GC Marker", "GC Marker", "GC Marker", "GC Marker", "GC Marker", 
                "GC Marker", "GC Marker", "GC Marker", 
                "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", 
                "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", 
                "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", 
                "Epithelial Marker", "Epithelial Marker", "Epithelial Marker", "Epithelial Marker")

heatmap_markers <- Heatmap(logcpm_markers_mat, cluster_rows = TRUE, cluster_columns = TRUE, 
                       name = "log2CPM + c", 
                       row_km = 3, 
                       rect_gp = gpar(col = "white", lwd = 0.25), 
                       clustering_distance_rows = "euclidean", 
                       row_names_gp = gpar(fontsize = 14), 
                       heatmap_legend_param = list(title_gp = gpar(fontsize = 14), 
                                                   labels_gp = gpar(fontsize = 14)), 
                       top_annotation = HeatmapAnnotation(df = data.frame(Category = sample_groups)))

