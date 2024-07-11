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

## load count file
count_file <- read.csv('E:/paper-files/novaseq_small_mevs_counts_final.csv', sep=',', header = TRUE)
colnames(count_file)[1] = "mirna" #change column 1 name because imports wonky
count_file <- column_to_rownames(count_file, "mirna") #changes the first column to rownames and now the df should match the metadata


#create vectors and dataframe containing metadata for the samples, however you can also load a metadatafile created in excel using the read.csv f(x) as well.
genotype <- c("ko", "ko", "ko", "ko", "ko", "ko", "ko", "ko", "ko", "wt", "wt", "wt", "wt", "wt",
              "wt", "wt", "wt", "wt", "wt", "wt", "ko", "ko", "ko", "ko", "ko", "ko", "ko", "ko",
              "ko", "ko", "wt", "wt", "wt", "wt", "wt", "wt", "wt", "wt", "wt")
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

#prefiltering - recommended but not required, creates smaller object to increase speed of computation
#this filters out genes that have greater or equal to 10 counts in at least '7' samples or whatever smallest_group is set as
keep <- rowSums(counts(count_file_dds) >= 10) >= smallest_group 
count_file_dds <- count_file_dds[keep,]

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
