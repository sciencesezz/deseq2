######################################Deseq2 Human EVs

##change path that packages are loaded to
.libPaths("E:/R-packages") #make sure this is set up locally
.libPaths()#prints out the paths, in position 1 is where the packages will be loaded

options(lib.loc = "E:/R-packages")

package_load <- c("tidyverse", "DESeq2", "tibble", "dplyr", "tidyr", "readr", "stringr", 
                  "ggplot2", "tidybulk", "ComplexHeatmap", "tidyHeatmap", "ggrepel", "plotly", 
                  "RColorBrewer", "pheatmap", "apeglm", "ashr", "annotables", "edgeR")
lapply(package_load, require, character.only = TRUE)

##create required variables
smallest_group <- 3 #smallest group for pre-filtering, change depending on dataset
FC <- 0.32 # Fold-change cutoff DESeq analysis
FDR <- 0.05 # FDR cutoff for DESeq analysis
alpha <- 0.1 # independent filtering, default for DESeq analysis

## load count file
count_file <- read.csv('E:/paper-files/novaseq_small_hevs_counts_xhgsocxsbt.csv', sep=',', header = TRUE)
colnames(count_file)[1] = "mirna" #change column 1 name because imports wonky
count_file <- column_to_rownames(count_file, "mirna") #changes the first column to rownames and now the df should match the metadata


#create vectors and dataframe containing metadata for the samples, however you can also load a metadatafile created in excel using the read.csv f(x) as well.
condition <- c("benign", "benign", "benign", "benign", "benign", "HGSOC", "HGSOC", "HGSOC", 
               "normal", "normal", "normal", "normal", "SBT", "SBT", "SBT", "SBT", "SBT")
group <- c("control", "control", "control", "control", "control", "HGSOC", "HGSOC", "HGSOC", 
           "control", "control", "control", "control", "SBT", "SBT", "SBT", "SBT", "SBT")
metadata <- data.frame(condition, group)

rownames(metadata) <- c("Benign_rep1", "Benign_rep2", "Benign_rep3", "Benign_rep4", "Benign_rep5",  
                        "HGSOC_rep1", "HGSOC_rep3", "HGSOC_rep4", "Normal_rep1", "Normal_rep2", "Normal_rep3", "Normal_rep4", 
                        "SBT_rep1",  "SBT_rep3",  "SBT_rep4",  "SBT_rep5",  "SBT_rep6")


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

write.csv(count_file_cpm, "E:/paper-files/hevs_logcpmc.csv", row.names = TRUE)


## Creation of the DESeqDataSet to include metadata information
count_file_dds <- DESeq2::DESeqDataSetFromMatrix(countData = count_file,
                                                 colData = metadata,
                                                 design = ~group)
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

write.csv(results_sig, "E:/HGSOC_v_control_hevs_mirna_deseq_xlfcshrink.csv", row.names = FALSE)
write.csv(results2_sig, "E:/SBT_v_control_hevs_mirna_deseq_xlfcshrink.csv", row.names = FALSE)
