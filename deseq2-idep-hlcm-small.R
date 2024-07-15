######################################Deseq2 Human Small LCM############################

##change path that packages are loaded to
.libPaths("E:/R-packages") #make sure this is set up locally
.libPaths()#prints out the paths, in position 1 is where the packages will be loaded

options(lib.loc = "E:/R-packages")

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')

if (!requireNamespace('circlize', quietly = TRUE))
  BiocManager::install('circlize')
require(circlize)

if (!requireNamespace('digest', quietly = TRUE))
  BiocManager::install('digest')
require(digest)

if (!requireNamespace('cluster', quietly = TRUE))
  BiocManager::install('cluster')
require(cluster)

package_load <- c("tidyverse", "DESeq2", "tibble", "dplyr", "tidyr", "readr", "stringr", 
                  "ggplot2", "tidybulk", "ComplexHeatmap", "tidyHeatmap", "ggrepel", "plotly", 
                  "RColorBrewer", "pheatmap", "apeglm", "ashr", "annotables", "edgeR")
lapply(package_load, require, character.only = TRUE)

##create required variables
smallest_group <- 3 #smallest group for pre-filtering, change depending on dataset
FC <- 0.32 # Fold-change cutoff DESeq analysis
FDR <- 0.05 # FDR cutoff for DESeq analysis
alpha <- 0.1 # independent filtering, default for DESeq analysis
condition_colours <- c("HGSOC" = "#D81B60", "Normal" = "#1E88E5", "SBT" = "#004D40", "Benign" = "#5D286B")
group_colours <- c("Control" = "#1E88E5", "SBT" = "#004D40", "HGSOC" = "#D81B60")

## load count file
count_file <- read.csv('E:/paper-files/novaseq_small_hlcm_counts_final_hgsocx.csv', sep=',', header = TRUE)
colnames(count_file)[1] = "mirna" #change column 1 name because imports wonky
count_file <- column_to_rownames(count_file, "mirna") #changes the first column to rownames and now the df should match the metadata


#create vectors and dataframe containing metadata for the samples, however you can also load a metadatafile created in excel using the read.csv f(x) as well.
condition <- c("Benign", "Benign", "Benign", "Benign", "Benign", "HGSOC", "HGSOC", "HGSOC", 
               "Normal", "Normal", "Normal", "Normal", "SBT", "SBT", "SBT", "SBT", "SBT", "SBT")
group <- c("Control", "Control", "Control", "Control", "Control", "HGSOC", "HGSOC", "HGSOC", 
           "Control", "Control", "Control", "Control", "SBT", "SBT", "SBT", "SBT", "SBT", "SBT")
metadata <- data.frame(condition, group)

rownames(metadata) <- c("Benign_rep1", "Benign_rep2", "Benign_rep3", "Benign_rep4", "Benign_rep5",  
                        "HGSOC_rep1", "HGSOC_rep3", "HGSOC_rep4", "Normal_rep1", "Normal_rep2", "Normal_rep3", "Normal_rep4", 
                        "SBT_rep1", "SBT_rep2", "SBT_rep3",  "SBT_rep4",  "SBT_rep5",  "SBT_rep6")


#check metadata file and count file match up
all(rownames(metadata) == colnames(count_file))

#Filter data, at least a count of 10 in the smallest sample group specified above
keep <- rowSums(count_file >= 10) >= smallest_group 
count_file <- count_file[keep,]

#convert to cpm to normalise for library size using edgeR so you have logCPM + c values to graph
count_file_cpm <- edgeR::cpm(count_file)
print(count_file_cpm)

#log2 transformation of the cpm + c (c is the pseudo count so that log2 isn't performed on 0, so 0 = 2)
pseudo <- 4
count_file_cpm <- log2(count_file_cpm + pseudo)
rownames(count_file_cpm) <- rownames(count_file)


write.csv(count_file_cpm, "E:/paper-files/hlcm_small_logcpmc.csv", row.names = TRUE)


## Creation of the DESeqDataSet to include metadata information
count_file_dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_file,
                                                 colData = metadata,
                                                 design = ~group)

#set the factor level, tell Deseq which level to compare against
count_file_dds$group <- relevel(count_file_dds$group, ref = "Control")
order_condition <- c("Normal", "Benign", "SBT", "HGSOC")
order_group <- c("Control", "SBT", "HGSOC")
count_file_dds$condition <- factor(count_file_dds$condition, levels = order_condition)
count_file_dds$group <- factor(count_file_dds$group, levels = order_group)
########here VST transformation for PCA and corr plot

##data transformation for PCA and Corr Plots of data
#you can choose multiple options, but i will go with VST in the deseq package
count_file_dds_vst <- varianceStabilizingTransformation(count_file_dds, blind = TRUE)

#so generate the PCA plot manually with ggplot
pcaData <- plotPCA(count_file_dds_vst, intgroup = c("condition"), 
                   returnData = TRUE)

percentVar <- round(100 * attr(pcaData, "percentVar"))

pcaplot <- ggplot(pcaData, aes(x = PC1, y = PC2, colour = condition)) +
  geom_point(size = 3, alpha = 0.7) + 
  scale_color_manual(values = condition_colours) +  
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance"))

print(pcaplot)

ggsave(filename = "E:/paper-files/hlcm_small_pca_small.png", plot = pcaplot, width = 4, height = 3, dpi = 800)
ggsave(filename = "E:/paper-files/hlcm_small_pca_big.png", plot = pcaplot, width = 8, height = 6, dpi = 800)

#plot correlation (non complex heatmap)
count_file_mat_vst <- assay(count_file_dds_vst) #extract the vst matrix from the object
corr_value <- cor(count_file_mat_vst) #compute pairwise correlation values

#######################ComplexHeatmap Correlation Plot############################
####making a heatmap using the ComplexHeatmap function
#need to make the annotation bars for heatmap annotation
#top annotation
#relevel the metadata to match the PCA plot
metadata$condition <- factor(metadata$condition, levels = order_condition)
metadata$group <- factor(metadata$group, levels = order_group)

ha_top <- HeatmapAnnotation(
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
ha_row <- HeatmapAnnotation(which = "row",
                            condition = metadata$condition,
                            col = list(condition = condition_colours 
                            ),
                            annotation_height = 0.3, 
                            annotation_width = unit(1, 'cm'), 
                            gap = unit(1, 'mm'), 
                            border = TRUE, 
                            show_legend = FALSE, 
                            show_annotation_name = FALSE, 
                            annotation_legend_param = list(
                              condition = list(
                                nrow = 4, 
                                title = "Condition", 
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

##converted counts similar to iDEP
count_file_dds <- DESeq2::estimateSizeFactors(count_file_dds)
sizeFactors(count_file_dds)

dds <- DESeq2::DESeq(count_file_dds)

results <- DESeq2::results(dds, 
                           contrast = c("group", "HGSOC", "control"), #name of the factor, name of the numerator, then denominator
                           alpha = alpha,
                           independentFiltering = TRUE) 

results2 <- DESeq2::results(dds, 
                            contrast = c("group", "SBT", "control"), #name of the factor, name of the numerator, then denominator
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

write.csv(results_sig, "E:/paper-files/HGSOC_v_control_hlcm_small_deseq_xlfcshrink.csv", row.names = FALSE)
write.csv(results2_sig, "E:/paper-files/SBT_v_control_hlcm_small_deseq_xlfcshrink.csv", row.names = FALSE)
