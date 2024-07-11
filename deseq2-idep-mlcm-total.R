######################################Deseq2 Human Total LCM############################

##change path that packages are loaded to
.libPaths("E:/R-packages") #make sure this is set up locally
.libPaths()#prints out the paths, in position 1 is where the packages will be loaded

options(lib.loc = "E:/R-packages")

package_load <- c("tidyverse", "DESeq2", "tibble", "dplyr", "tidyr", "readr", "stringr", 
                  "ggplot2", "tidybulk", "ComplexHeatmap", "tidyHeatmap", "ggrepel", "plotly", 
                  "RColorBrewer", "pheatmap", "apeglm", "ashr", "annotables", "edgeR")
lapply(package_load, require, character.only = TRUE)

##create required variables
smallest_group <- 7 #smallest group for pre-filtering, change depending on dataset
FC <- 0.32 # Fold-change cutoff DESeq analysis
FDR <- 0.05 # FDR cutoff for DESeq analysis
alpha <- 0.1 # independent filtering, default for DESeq analysis
## load count file
count_file <- read.csv('E:/paper-files/novaseq_total_mlcm_all_counts_final.csv', sep=',', header = TRUE)
colnames(count_file)[1] = "ensgene" #change column 1 name because imports wonky
count_file <- column_to_rownames(count_file, "ensgene") #changes the first column to rownames and now the df should match the metadata

#create vectors and dataframe containing metadata for the samples, however you can also load a metadatafile created in excel using the read.csv f(x) as well.
genotype <- c("ko", "ko", "ko", "ko", "ko", "ko", "ko", "ko", "ko", "ko", "ko", "ko", "ko", "ko",
              "ko", "ko", "ko", "ko", "ko", "ko", "ko", "wt", "wt", "wt", "wt", "wt", "wt", "wt",
              "wt", "wt", "wt", "wt", "wt", "wt", "wt", "wt", "wt", "wt", "wt", "wt", "wt", "wt")
condition <- c("adenoma", "adenoma", "adenoma", "adenoma", "adenoma", "adenoma", "adenoma",
               "stroma", "stroma", "stroma", "stroma", "stroma", "stroma", "stroma",
               "sex cords", "sex cords", "sex cords", "sex cords", "sex cords", "sex cords", "sex cords",
               "mature GCs", "mature GCs", "mature GCs", "mature GCs", "mature GCs",  "mature GCs",  "mature GCs", 
               "primary GCs", "primary GCs", "primary GCs", "primary GCs", "primary GCs", "primary GCs", "primary GCs", 
               "stroma", "stroma", "stroma", "stroma", "stroma", "stroma", "stroma")
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
##data transformation for PCA and Corr Plots of data
#you can choose multiple options, but i will go with VST in the deseq package


##converted counts similar to iDEP
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

