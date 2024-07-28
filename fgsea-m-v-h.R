rm(list = ls())
##change path that packages are loaded to
.libPaths("E:/R-packages") #make sure this is set up locally
.libPaths()#prints out the paths, in position 1 is where the packages will be loaded

options(lib.loc = "E:/R-packages")

package_load <- c("tidyverse", "DESeq2", "tibble", "dplyr", "tidyr", "readr", "stringr", 
                  "ggplot2", "tidybulk", "ComplexHeatmap", "tidyHeatmap", "ggrepel", "plotly", 
                  "RColorBrewer", "pheatmap", "apeglm", "ashr", "annotables", "edgeR", "VennDiagram", "limma", 
                  "tools")
lapply(package_load, require, character.only = TRUE)
library(clusterProfiler)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(enrichplot)
library(stringr)

#############################GSEA########################################

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("fgsea")
install.packages("data.table")
library(data.table)
library(fgsea)

#prepare variables
min_size <- 15
max_size <- 500
padj_cutoff <- 0.05

#prepare functions
#prepare gmt
#this function allows to filter gmt files by a background gene list (filtered gene list pre dedeq2)

## Function: Adjacency matrix to list -------------------------
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}
## Function: prepare_gmt --------------------------------------
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE){
  # for debug
  #file <- gmt_files[1]
  #genes_in_data <- df$gene_symbol
  
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with the genes as rows and for each go annotation (columns) the values are 0 or 1
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  #Subset to the genes that are present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]] # filter for gene sets with more than 5 genes annotated
  # And get the list again
  final_list <- matrix_to_list(mat) # for this we use the function we previously defined
  
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('Wohoo! .gmt conversion successfull!:)')
  return(final_list)
}

mygenes_mouse <- read.csv('E:/paper-files/mlcm_total_logcpmc.csv', sep=',', header = TRUE)
colnames(mygenes_mouse)[1] <- "ensgene"
mygenes_mouse <- left_join(x = mygenes_mouse,
                           y = grcm38 [, c("ensgene", "symbol", "entrez", "biotype", "description")], 
                           by = "ensgene")
mygenes_mouse <- subset(mygenes_mouse, biotype == "protein_coding")
mygenes_mouse <- mygenes_mouse$symbol

mygenes_humans <- read.csv('E:/paper-files/hlcm_total_logcpmc.csv', sep=',', header = TRUE)
colnames(mygenes_humans)[1] <- "ensgene"
mygenes_humans <- left_join(x = mygenes_humans,
                            y = grch38 [, c("ensgene", "symbol", "entrez", "biotype", "description")], 
                            by = "ensgene")
#get results
write.csv(mygenes_humans, "E:/paper-files/hlcm_total_logcpmc_anno.csv", row.names = FALSE)

#set up filtered gmt pathways file
go_filtered_gmt <- prepare_gmt("E:/paper-files/m5.go.v2023.2.Mm.symbols.gmt", 
                               mygenes_mouse, 
                               savefile = TRUE)

#get ranked gene list
sxcd_go <- as.data.frame(sxcd_lfc)
sxcd_go <- sxcd_go %>%
  mutate(symbol = str_to_title(symbol))
sxcd_go$log2FoldChange <- as.numeric(sxcd_go$log2FoldChange)
#sxcd_go$padj <- as.numeric(sxcd_go$padj)
sxcd_go <- sxcd_go%>%
  arrange(dplyr::desc(log2FoldChange))
sxcd_go <- sxcd_go[,-3]
sxcd_go <- setNames(sxcd_go$log2FoldChange, sxcd_go$symbol)

#alternative way of ranking
#rankings <- sign(sxcd_go$log2FoldChange)*(-log10(sxcd_go$padj))
#names(rankings) <- sxcd_go$symbol

#alternative way of importing gmt files
#go_pathway <- gmtPathways("E:/paper-files/m5.go.v2023.2.Mm.symbols.gmt")
#tumour_pathway <- gmtPathways("E:/paper-files/m5.mpt.v2023.2.Mm.symbols.gmt")
#cc_pathway <- gmtPathways("E:/paper-files/m5.go.cc.v2023.2.Mm.symbols.gmt")


# Run fgsea
fgsea_sxcd_go <- fgsea(pathways = go_filtered_gmt, 
                       stats = sxcd_go ,
                       minSize = min_size,
                       maxSize = max_size, 
                       scoreType = 'std') #both positive and negative rankings

fgsea_sxcd_cc <- fgseaMultilevel(pathways = cc_pathway, 
                                 stats = rankings ,
                                 minSize = min_size,
                                 maxSize = max_size)

fgsea_sxcd_tumour <- fgsea(pathways = tumour_pathway, 
                           stats = sxcd_go ,
                           minSize = min_size,
                           maxSize = max_size)

