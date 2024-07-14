######################################Deseq2 Human Total LCM############################

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
smallest_group <- 7 #smallest group for pre-filtering, change depending on dataset
FC <- 0.32 # Fold-change cutoff DESeq analysis
FDR <- 0.05 # FDR cutoff for DESeq analysis
alpha <- 0.1 # independent filtering, default for DESeq analysis
condition_colours <- c("Adenoma" = "indianred", "Mature GCs" = "hotpink", 
                       "Primary GCs" = "springgreen3", "Sex cords" = "steelblue2", 
                       "Stroma" = "goldenrod")
genotype_colours <- c("KO" = "black", "WT" = "grey57")

## load count file
count_file <- read.csv('E:/paper-files/novaseq_total_mlcm_all_counts_final.csv', sep=',', header = TRUE)
colnames(count_file)[1] = "ensgene" #change column 1 name because imports wonky
count_file <- column_to_rownames(count_file, "ensgene") #changes the first column to rownames and now the df should match the metadata

#create vectors and dataframe containing metadata for the samples, however you can also load a metadatafile created in excel using the read.csv f(x) as well.
genotype <- c("KO", "KO", "KO", "KO", "KO", "KO", "KO", "KO", "KO", "KO", "KO", "KO", "KO", "KO",
              "KO", "KO", "KO", "KO", "KO", "KO", "KO", "WT", "WT", "WT", "WT", "WT", "WT", "WT",
              "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT", "WT")
condition <- c("Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma",
               "Stroma", "Stroma", "Stroma", "Stroma", "Stroma", "Stroma", "Stroma",
               "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords",
               "Mature GCs", "Mature GCs", "Mature GCs", "Mature GCs", "Mature GCs",  "Mature GCs",  "Mature GCs", 
               "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs", 
               "Stroma", "Stroma", "Stroma", "Stroma", "Stroma", "Stroma", "Stroma")
group <- c("Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma", 
           "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
           "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", 
           "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
           "Control", "Control", "Control", "Control", "Control", "Control", "Control", 
           "Control", "Control", "Control", "Control", "Control", "Control", "Control")
input <- c("high", "high", "high", "high", "high", "high", "high", 
           "high", "high", "high", "high", "high", "high", "high", 
           "low", "low", "low", "low", "low", "low", "low", 
           "high", "high", "high", "high", "high", "high", "high", 
           "low", "low", "low", "low", "low", "low", "low", 
           "high", "high", "high", "high", "high", "high", "high")
metadata <- data.frame(genotype, condition, group, input)
rownames(metadata) <- c("Adenoma_1", "Adenoma_2", "Adenoma_3", "Adenoma_4", "Adenoma_5", "Adenoma_6", "Adenoma_7", 
                        "KO_Stro_1", "KO_Stro_2", "KO_Stro_3", "KO_Stro_4", "KO_Stro_5", "KO_Stro_6", "KO_Stro_7", 
                        "Sxcd_1", "Sxcd_2", "Sxcd_3", "Sxcd_4", "Sxcd_5", "Sxcd_6", "Sxcd_7", 
                        "GC_1", "GC_2", "GC_3", "GC_4", "GC_5", "GC_6", "GC_7", 
                        "Prim_1", "Prim_2", "Prim_3", "Prim_4", "Prim_5", "Prim_6", "Prim_7", 
                        "WT_Stro_1", "WT_Stro_2", "WT_Stro_3", "WT_Stro_4", "WT_Stro_5", "WT_Stro_6", "WT_Stro_7")

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
write.csv(count_file_cpm, "E:/paper-files/mlcm_total_logcpmc.csv", row.names = TRUE)

## Creation of the DESeqDataSet to include metadata information
count_file_dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_file,
                                                 colData = metadata,
                                                 design = ~group)

#set the factor level, tell Deseq which level to compare against
count_file_dds$group <- relevel(count_file_dds$group, ref = "Control")
order_condition <- c("Adenoma", "Sex cords", "Primary GCs", "Mature GCs", "Stroma")
order_geno <- c("WT", "KO")
count_file_dds$condition <- factor(count_file_dds$condition, levels = order_condition)
count_file_dds$genotype <- factor(count_file_dds$genotype, levels = order_geno)

##data transformation for PCA and Corr Plots of data
#you can choose multiple options, but i will go with VST in the deseq package
count_file_dds_vst <- vst(count_file_dds, blind = TRUE)

#plot PCA with Deseq2, but you can't do two groups
DESeq2::plotPCA(count_file_dds_vst, intgroup = c("condition"))

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

ggsave(filename = "E:/paper-files/mlcm_total_pca_small.png", plot = pcaplot, width = 4, height = 3, dpi = 800)
ggsave(filename = "E:/paper-files/mlcm_total_pca_big.png", plot = pcaplot, width = 8, height = 6, dpi = 800)


#plot correlation (non complex heatmap)
count_file_mat_vst <- assay(count_file_dds_vst) #extract the vst matrix from the object
corr_value <- cor(count_file_mat_vst) #compute pairwise correlation values

corrplot <- pheatmap(corr_value, annotation = select(metadata, condition), 
         fontsize_row = 8, 
         fontsize_col = 8, 
         main = "Correlation Values of Sample Transcriptomes")
print(corrplot)
ggsave(filename = "E:/paper-files/mlcm_total_corr_small.png", plot = corrplot, width = 4, height = 3, dpi = 800)
ggsave(filename = "E:/paper-files/mlcm_total_corr_big.png", plot = corrplot, width = 8, height = 6, dpi = 800)

#######################ComplexHeatmap Correlation Plot############################
####making a heatmap using the ComplexHeatmap function
#need to make the annotation bars for heatmap annotation
#top annotation
#relevel the metadata to match the PCA plot
metadata$condition <- factor(metadata$condition, levels = order_condition)
metadata$genotype <- factor(metadata$genotype, levels = order_geno)

ha_top <- HeatmapAnnotation(
  condition = metadata$condition,
  genotype = metadata$genotype,
  col = list(condition = condition_colours, 
             genotype = genotype_colours 
  ),
  annotation_height = unit(3, "mm"), 
  annotation_width = unit(0.5, 'cm'), 
  gap = unit(0.5, 'mm'), 
  border = TRUE, 
  annotation_legend_param = list(
    condition = list(
      nrow = 5, 
      title = "Condition", 
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
                            condition = metadata$condition,
                            col = list(genotype = genotype_colours, 
                                       condition = condition_colours 
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
                              condition = list(
                                nrow = 5, 
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

########################################Deseq2 DEGS###############################
##then do the actual differential gene expression analysis
count_file_dds <- DESeq2::estimateSizeFactors(count_file_dds)
sizeFactors(count_file_dds)
dds <- DESeq2::DESeq(count_file_dds)


results <- DESeq2::results(dds, 
                           contrast = c("group", "Adenoma", "Control"), #name of the factor, name of the numerator, then denominator
                           alpha = alpha,
                           independentFiltering = TRUE) 

results2 <- DESeq2::results(dds, 
                            contrast = c("group", "Sex cords", "Control"), #name of the factor, name of the numerator, then denominator
                            alpha = alpha, 
                            independentFiltering = TRUE)

#save the significant results
results_all <- data.frame(results)
results_sig <- subset(results_all, padj < 0.05)
results2_all <- data.frame(results2)
results2_sig <- subset(results2_all, padj < 0.05)

## annotate significant genes
results_anno <- data.frame(results_sig) 
results_anno <- rownames_to_column(results_anno, var = "ensgene")
results_anno <- left_join(x = results_anno,
                          y = grcm38 [, c("ensgene", "symbol", "entrez", "biotype", "description")], 
                          by = "ensgene")
results2_anno <- data.frame(results2_sig)
results2_anno <- rownames_to_column(results2_anno, var = "ensgene")
results2_anno <- left_join(x = results2_anno,
                           y = grcm38 [, c("ensgene", "symbol", "entrez", "biotype", "description")], 
                           by = "ensgene")

#get results
write.csv(results_anno, "E:/paper-files/mlcm_total_adeno_v_control_deseq_xlfcshrink.csv", row.names = FALSE)
write.csv(results2_anno, "E:/paper-files/mlcm_total_sxcd_v_control_deseq_xlfcshrink.csv", row.names = FALSE)

