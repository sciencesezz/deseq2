##############################################DESeq2 Mouse######################a

##change path that packages are loaded to
.libPaths("E:/R-packages") #make sure this is set up locally
.libPaths()#prints out the paths, in position 1 is where the packages will be loaded

options(lib.loc = "E:/R-packages")

package_load <- c("tidyverse", "DESeq2", "tibble", "dplyr", "tidyr", "readr", "stringr", 
                  "ggplot2", "tidybulk", "ComplexHeatmap", "tidyHeatmap", "ggrepel", "plotly", 
                  "RColorBrewer", "pheatmap", "apeglm", "ashr", "annotables", "edgeR")
lapply(package_load, require, character.only = TRUE)


##create required variables
smallest_group <- 9 #smallest group for pre-filtering, change depending on dataset
FC <- 0.32 # Fold-change cutoff DESeq analysis
FDR <- 0.1 # FDR cutoff for DESeq analysis
alpha <- 0.1 # independent filtering, default for DESeq analysis
group_colours <- c("Adenoma" = "#D81B60", 
                   "Control" = "#1E88E5", 
                   "Sex cords" = "#004D40")     #colour blind safe
genotype_colours <- c("KO" = "black", "WT" = "grey57")
age_colours <- c("1_year" = "#D81B60", "3_months" = "#004D40")
condition_colours <- c("1y_KO" = "#D81B60", "1y_WT" = "#1E88E5", "3m_KO" = "#004D40", "3m_WT" = "turquoise")

## load count file
count_file <- read.csv('E:/paper-files/novaseq_small_mevs_counts_final.csv', sep=',', header = TRUE)
colnames(count_file)[1] = "mirna" #change column 1 name because imports wonky
count_file <- column_to_rownames(count_file, "mirna") #changes the first column to rownames and now the df should match the metadata


#create vectors and dataframe containing metadata for the samples, however you can also load a metadatafile created in excel using the read.csv f(x) as well.
genotype <- c("KO", "KO", "KO", "KO", "KO", "KO", "KO", "KO", "KO", "WT", "WT", "WT", "WT", "WT",
              "WT", "WT", "WT", "WT", "WT", "WT", "KO", "KO", "KO", "KO", "KO", "KO", "KO", "KO",
              "KO", "KO", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT")
age <- c("1_year", "1_year", "1_year", "1_year", "1_year", "1_year", "1_year",
         "1_year", "1_year", "1_year", "1_year", "1_year", "1_year", "1_year",
         "1_year", "1_year", "1_year", "1_year", "1_year", "1_year", "3_months",
         "3_months", "3_months", "3_months", "3_months", "3_months",  "3_months",  "3_months", 
         "3_months", "3_months", "3_months", "3_months", "3_months", "3_months", "3_months", 
         "3_months", "3_months", "3_months", "3_months")
group <- c("Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma", 
           "Adenoma", "Adenoma", "Control", "Control", "Control", "Control", "Control", 
           "Control", "Control", "Control", "Control", "Control", "Control", "Sex cords", 
           "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", 
           "Sex cords", "Sex cords", "Control", "Control", "Control", "Control", "Control", 
           "Control", "Control", "Control", "Control")
condition <- c("1y_KO", "1y_KO", "1y_KO", "1y_KO", "1y_KO", "1y_KO", "1y_KO", 
               "1y_KO", "1y_KO", "1y_WT", "1y_WT", "1y_WT", "1y_WT", "1y_WT", 
               "1y_WT", "1y_WT", "1y_WT", "1y_WT", "1y_WT", "1y_WT", "3m_KO", 
               "3m_KO", "3m_KO", "3m_KO", "3m_KO", "3m_KO", "3m_KO", "3m_KO", 
               "3m_KO", "3m_KO", "3m_WT", "3m_WT", "3m_WT", "3m_WT", "3m_WT", 
               "3m_WT", "3m_WT", "3m_WT", "3m_WT")
metadata <- data.frame(genotype, age, group, condition)

rownames(metadata) <- c("KO_1y_1", "KO_1y_2", "KO_1y_3",	"KO_1y_4",	"KO_1y_5",	"KO_1y_6",	"KO_1y_7",	"KO_1y_8",	"KO_1y_9",
                        "WT_1y_1",	"WT_1y_2", "WT_1y_3",	"WT_1y_4",	"WT_1y_5",	"WT_1y_6",	"WT_1y_7", "WT_1y_8",	"WT_1y_9",	"WT_1y_10",	"WT_1y_11",	
                        "KO_3m_1",	"KO_3m_2",	"KO_3m_3",	"KO_3m_4",	"KO_3m_5",	"KO_3m_6",	"KO_3m_7",	"KO_3m_8",	"KO_3m_9",	"KO_3m_10",	
                        "WT_3m_1",	"WT_3m_2",	"WT_3m_3",	"WT_3m_4",	"WT_3m_5",	"WT_3m_6",	"WT_3m_7",	"WT_3m_8",	"WT_3m_9")

#check metadata file and count file match up
all(rownames(metadata) == colnames(count_file))

#log2 transformation of the data
#Filter data, at least a count of 10 in the smallest sample group specified above
keep <- rowSums(count_file >= 10) >= smallest_group 
count_file <- count_file[keep,]

#convert to cpm to normalise for library size using edgeR
count_file_cpm <- edgeR::cpm(count_file)

#log2 transformation of the cpm + c (c is the pseudo count so that log2 isn't performed on 0, so 0 = 2)
pseudo <- 4
count_file_cpm <- log2(count_file_cpm + pseudo)

write.csv(count_file_cpm, "E:/paper-files/mevs_logcpmc.csv", row.names = TRUE)

## Creation of the DESeqDataSet to include metadata information
count_file_dds <- DESeqDataSetFromMatrix(countData = count_file,
                                         colData = metadata,
                                         design = ~ condition)
#set the factor level, tell Deseq which level to compare against
#count_file_dds$group <- relevel(count_file_dds$group, ref = "Control")
order_condition <- c("3m_WT", "3m_KO", "1y_WT", "1y_KO")
order_geno <- c("WT", "KO")
order_group <- c("Control", "Adenoma", "Sex cords")
order_age <- c("3_months", "1_year")

count_file_dds$condition <- factor(count_file_dds$condition, levels = order_condition)
count_file_dds$genotype <- factor(count_file_dds$genotype, levels = order_geno)
count_file_dds$group <- factor(count_file_dds$group, levels = order_group)
count_file_dds$age <- factor(count_file_dds$age, levels = order_age)

##data transformation for PCA and Corr Plots of data
#you can choose multiple options, but i will go with VST in the deseq package
count_file_dds_vst <- varianceStabilizingTransformation(count_file_dds, blind = TRUE)
##############################PCA Plot############################################
#plot PCA with Deseq2, but you can't do two groups
#DESeq2::plotPCA(count_file_dds_vst, intgroup = c("condition"))

#so generate the PCA plot manually with ggplot
pcaData <- plotPCA(count_file_dds_vst, intgroup = c("condition", "genotype"), 
                   returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaplot <- ggplot(pcaData, aes(x = PC1, y = PC2, colour = condition, shape = genotype)) +
  geom_point(size = 3, alpha = 0.7) + 
  scale_color_manual(values = condition_colours) +  
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance"))

print(pcaplot)

ggsave(filename = "E:/paper-files/mevs_small_pca_big.png", plot = pcaplot, width = 8, height = 6, dpi = 800)
ggsave(filename = "E:/paper-files/mevs_small_pca_small.png", plot = pcaplot, width = 4, height = 3, dpi = 800)


#######################ComplexHeatmap Correlation Plot############################
#plot correlation using pheatmap from deseq2 package I think
count_file_mat_vst <- assay(count_file_dds_vst) #extract the vst matrix from the object
corr_value <- cor(count_file_mat_vst) #compute pairwise correlation values

####making a heatmap using the ComplexHeatmap function
#need to make the annotation bars for heatmap annotation
#top annotation
#relevel the metadata to match the PCA plot
metadata$condition <- factor(metadata$condition, levels = order_condition)
metadata$genotype <- factor(metadata$genotype, levels = order_geno)
metadata$age <- factor(metadata$age, levels = order_age)

ha_top <- HeatmapAnnotation(
  age = metadata$age,
  genotype = metadata$genotype,
  col = list(age = age_colours, 
             genotype = genotype_colours 
  ),
  annotation_height = unit(3, "mm"), 
  annotation_width = unit(0.5, 'cm'), 
  gap = unit(0.5, 'mm'), 
  border = TRUE, 
  annotation_legend_param = list(
    age = list(
      nrow = 2, 
      title = "Age", 
      title_position = 'topleft',
      legend_direction = 'horizontal',
      title_gp = gpar(fontsize =12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'plain')),
    genotype = list(
      nrow = 2,
      title = 'Genotype',
      title_position = 'topleft',
      legend_direction = 'horizontal',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'plain'))))

#row annotation (order of variables is th eorder of the annotation bars)
ha_row <- HeatmapAnnotation(which = "row",
                            genotype = metadata$genotype,                          
                            age = metadata$age,
                            col = list(genotype = genotype_colours, 
                                       age = age_colours 
                            ),
                            annotation_height = 0.3, 
                            annotation_width = unit(1, 'cm'), 
                            gap = unit(1, 'mm'), 
                            border = TRUE, 
                            show_legend = FALSE, 
                            show_annotation_name = FALSE, 
                            annotation_legend_param = list(
                              genotype = list(
                                nrow = 2,
                                title = 'Genotype',
                                title_position = 'topcenter',
                                legend_direction = 'vertical',
                                title_gp = gpar(fontsize = 12, fontface = 'plain'),
                                labels_gp = gpar(fontsize = 12, fontface = 'plain')), 
                              age = list(
                                nrow = 2, 
                                title = "Age", 
                                title_position = 'topcenter',
                                legend_direction = 'vertical',
                                title_gp = gpar(fontsize =12, fontface = 'plain'),
                                labels_gp = gpar(fontsize = 12, fontface = 'plain')) 
                            ))



corrplot2 <- ComplexHeatmap::Heatmap(corr_value, 
                                     name = "Correlation of Expressed miRNAs", 
                                     top_annotation = ha_top,
                                     right_annotation = ha_row, 
                                     show_row_names = FALSE, 
                                     show_column_names = FALSE, 
                                     row_names_gp = gpar(fontsize = 10), 
                                     heatmap_legend_param = list(
                                       color_bar = "continuous", 
                                       legend_direction = "vertical", 
                                       title = "Correlation",
                                       title_gp = gpar(fontsize = 12, fontface = "bold"),
                                       labels_gp = gpar(fontsize = 12, fonface = "plain"),
                                       grid_width = unit(3, "mm"),
                                       grid_height = unit(3, "mm")), 
                                     rect_gp = gpar(col = "grey10", lwd = 0.5),
                                     cell_fun = function(j, i, x, y, width, height, fill) {
                                       grid.rect(x = x, y = y, width = width, height = height, 
                                                 gp = gpar(fill = NA, col = "grey10", lwd = 0.5))
                                     })


#save from viewer...
#png(file = "E:/paper-files/mlcm_small_corr.png", width = 1200 , height = 800, res = 800)

#draw the heatmap
draw(corrplot2,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE)

# Close the graphics device
#dev.off()
#######################ComplexHeatmap Log(cpm+c) data############################

ha_top2 <- HeatmapAnnotation(
  condition = metadata$condition,
  col = list(condition = condition_colours),
  annotation_height = unit(3, "mm"), 
  annotation_width = unit(0.5, 'cm'), 
  gap = unit(0.5, 'mm'), 
  border = TRUE, 
  annotation_legend_param = list(
    condition = list(
      nrow = 4, 
      title = "Condition", 
      title_position = 'topleft',
      legend_direction = 'horizontal',
      title_gp = gpar(fontsize =12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'plain'))))
  
#row annotation (order of variables is th eorder of the annotation bars)

corrplot3 <- ComplexHeatmap::Heatmap(count_file_cpm, 
                                     name = "Correlation of Expressed miRNAs", 
                                     top_annotation = ha_top2,
                                     show_row_names = FALSE, 
                                     show_column_names = FALSE, 
                                     row_names_gp = gpar(fontsize = 10),
                                     heatmap_legend_param = list(
                                       color_bar = "continuous", 
                                       legend_direction = "vertical", 
                                       title = "Log2(cpm + C)",
                                       title_gp = gpar(fontsize = 12, fontface = "bold"),
                                       labels_gp = gpar(fontsize = 12, fonface = "plain"),
                                       grid_width = unit(3, "mm"),
                                       grid_height = unit(3, "mm")))

#save from viewer...
#png(file = "E:/paper-files/mlcm_small_corr.png", width = 1200 , height = 800, res = 800)

#draw the heatmap
draw(corrplot3,
     heatmap_legend_side = "right",
     annotation_legend_side = "right",
     merge_legend = TRUE)


##converted counts similar to iDEP
count_file_dds <- DESeq2::estimateSizeFactors(count_file_dds)

dds <- DESeq(count_file_dds)

results <- results(dds, 
                   contrast = c("condition", "3m_KO", "3m_WT"), #name of the factor, name of the numerator, then denominator
                   alpha = alpha,
                   independentFiltering = TRUE)  #define FC variable at the beginning
results2 <- results(dds, 
                    contrast = c("condition", "1y_KO", "1y_WT"), #name of the factor, name of the numerator, then denominator
                    alpha = alpha,
                    independentFiltering = TRUE) 

results_all <- data.frame(results)
results_sig <- subset(results_all, padj < 0.05)
results_sig <- data.frame(results_sig)
results_sig <- rownames_to_column(results_sig, var = "mirna")

results2_all <- data.frame(results2)
results2_sig <- subset(results2_all, padj < 0.05)
results2_sig <- data.frame(results2_sig)
results2_sig <- rownames_to_column(results2_sig, var = "mirna")

#get results

write.csv(results_sig, "E:/paper-files/3m_KO_v_3m_WT_mevs_mirna_deseq_xlfc2.csv", row.names = FALSE)
write.csv(results2_sig, "E:/paper-files/1y_KO_v_1y_WT_mevs_mirna_deseq_xlfc2.csv", row.names = FALSE)
