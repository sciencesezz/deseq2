
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


#############################Upload files as objects##########################
sxcd <- read.csv('E:/paper-files/mlcm_total_sxcd_v_control_deseq_xlfcshrink.csv', sep=',', header = TRUE)
adeno <- read.csv('E:/paper-files/mlcm_total_adeno_v_control_deseq_xlfcshrink.csv', sep=',', header = TRUE)
sbt <- read.csv('E:/paper-files/hlcm_total_SBT_v_control_deseq_xlfcshrink.csv', sep=',', header = TRUE)
hgsoc <- read.csv('E:/paper-files/hlcm_total_HGSOC_v_control_deseq_xlfcshrink.csv', sep=',', header = TRUE)

#set all symbols into lower case so they can be matched, I've already done this, but for simplicity sake
#also keep this block of code instead of referring to casematch column as symbol.
sxcd$symbol <- tolower(sxcd$symbol)
adeno$symbol <- tolower(adeno$symbol)
sbt$symbol <- tolower(sbt$symbol)
hgsoc$symbol <- tolower(hgsoc$symbol)

#filter each dataset for mRNAs
sxcd_up <- filter(sxcd, log2FoldChange >= 0, biotype == "protein_coding", symbol != "")
sxcd_down <- filter(sxcd, log2FoldChange < 0, biotype == "protein_coding", symbol != "")
adeno_up <- filter(adeno, log2FoldChange >= 0, biotype == "protein_coding", symbol != "")
adeno_down <- filter(adeno, log2FoldChange < 0, biotype == "protein_coding", symbol != "")
sbt_up <- filter(sbt, log2FoldChange >= 0, biotype == "protein_coding", symbol != "")
sbt_down <- filter(sbt, log2FoldChange < 0, biotype == "protein_coding", symbol != "")
hgsoc_up <- filter(hgsoc, log2FoldChange >= 0, biotype == "protein_coding", symbol != "")
hgsoc_down <- filter(hgsoc, log2FoldChange < 0, biotype == "protein_coding", symbol != "")

##select just mRNA symbols for venn comparison
sxcd_up <- arrange(sxcd_up, -log2FoldChange) %>%
  dplyr::select(symbol)
sxcd_down <- arrange(sxcd_down, log2FoldChange) %>%
  dplyr::select(symbol)
adeno_up <- arrange(adeno_up, -log2FoldChange) %>%
  dplyr:: select(symbol)
adeno_down <- arrange(adeno_down, log2FoldChange) %>%
  dplyr::select(symbol)
sbt_up <- arrange(sbt_up, -log2FoldChange) %>%
  dplyr::select(symbol)
sbt_down <- arrange(sbt_down, log2FoldChange) %>%
  dplyr::select(symbol)
hgsoc_up <- arrange(hgsoc_up, -log2FoldChange) %>%
  dplyr::select(symbol)
hgsoc_down <- arrange(hgsoc_down, log2FoldChange) %>%
  dplyr::select(symbol)

#coerce those symbol lists to vector for venn analysis
sxcd_up <- sxcd_up$symbol
sxcd_down <- sxcd_down$symbol
adeno_up <- adeno_up$symbol
adeno_down <- adeno_down$symbol
sbt_up <- sbt_up$symbol
sbt_down <- sbt_down$symbol
hgsoc_up <- hgsoc_up$symbol
hgsoc_down <- hgsoc_down$symbol
##############################VENN##############################################
venn_data_ordered_up <- list(Sexcords = sxcd_up, 
                          SBT = sbt_up,
                          Adenoma = adeno_up, 
                          HGSOC = hgsoc_up)

venn_data_ordered_down <- list(Sexcords = sxcd_down, 
                             SBT = sbt_down,
                             Adenoma = adeno_down, 
                             HGSOC = hgsoc_down)

venn_up_mh <- venn.diagram(x = venn_data_ordered_up,
                             category.names = c("Sexcords", "SBT", "Adenoma", "HGSOC"), 
                             filename = "E:/paper-files/images/venn_up_mh_big.tiff",
                             disable.logging = FALSE, 
                             imagetype = "tiff", 
                             main = "Mouse versus Human", 
                             main.fontface = "bold", 
                             sub = "Upregulated DEGs", 
                             col = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                             fill = c("palegreen1", "mediumpurple1","hotpink", "grey"), 
                             alpha = 0.3, 
                          fontfamily = "sans", 
                          resolution = 800, 
                          width = 8, 
                          height = 8, 
                          units = "in", 
                          sub.fontface = "bold", 
                          cat.fontfamily = "sans", 
                          cat.fontface = "bold", 
                          cex = 3,
                          cat.cex = 2.15, 
                          bg = "transparent"
                          )


venn_down_mh <- venn.diagram(x = venn_data_ordered_down,
                           category.names = c("Sexcords", "SBT", "Adenoma", "HGSOC"), 
                           filename = "E:/paper-files/images/venn_down_mh_big.tiff",
                           disable.logging = FALSE, 
                           imagetype = "tiff", 
                           main = "Mouse versus Human", 
                           main.fontface = "bold", 
                           sub = "Downregulated DEGs", 
                           col = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                           fill = c("palegreen1", "mediumpurple1","hotpink", "grey"), 
                           alpha = 0.3, 
                           fontfamily = "sans", 
                           resolution = 800, 
                           width = 8, 
                           height = 8, 
                           units = "in", 
                           sub.fontface = "bold", 
                           cat.fontfamily = "sans", 
                           cat.fontface = "bold", 
                           cex = 3, 
                           cat.cex = 2.15)

##upset plots of common genes:
install.packages("UpSetR")
library(UpSetR)
install.packages("ggplotify")
library(ggplotify)
install.packages("Cairo")
library(Cairo)

CairoPNG("E:/paper-files/images/upset_genes_mh_up.png", width = 8, height = 8, units = "in", res = 800, 
         bg = "transparent")

upset(
  fromList(venn_data_ordered_up),
  order.by = "freq",  # Order by frequency
  point.size = 3.5,  # Size of points in the matrix
  line.size = 1,  # Size of lines connecting sets
  main.bar.color = "mediumpurple",  # Color of the main bar
  sets.bar.color = "gray",  # Color of the sets bar
  matrix.color = "red3",  # Color of the matrix points
  text.scale = c(1.8, 1.8, 1.2, 1.2, 1.8, 1.8),  # Scale of text elements
  sets.x.label = "No. of Genes in Set",  # Label for sets bar
  keep.order = TRUE,  # Keep the order of sets
  empty.intersections = "on"  # Show empty intersections
)
# Close the graphics device
dev.off()

CairoPNG("E:/paper-files/images/upset_genes_mh_down.png", width = 8, height = 8, units = "in", res = 800)

upset(
  fromList(venn_data_ordered_down),
  order.by = "freq",  # Order by frequency
  point.size = 3.5,  # Size of points in the matrix
  line.size = 1,  # Size of lines connecting sets
  main.bar.color = "darkgreen",  # Color of the main bar
  sets.bar.color = "gray",  # Color of the sets bar
  matrix.color = "royalblue4",  # Color of the matrix points
  text.scale = c(1.8, 1.8, 1.2, 1.2, 1.8, 1.8),  # Scale of text elements
  sets.x.label = "No. of Genes in Set",  # Label for sets bar
  keep.order = TRUE,  # Keep the order of sets
  empty.intersections = "on"  # Show empty intersections
)
# Close the graphics device
dev.off()

# Create a list of vectors
vectors_list_up <- list(sxcd_up, sbt_up, adeno_up, hgsoc_up)
elements_up <- calculate.overlap(vectors_list_up)
# Find the maximum length of the vectors
max_length_up <- max(sapply(elements_up, length))
# Pad the vectors with NA to make them of equal length
padded_list_up <- lapply(elements_up, function(vec) {
  length(vec) <- max_length_up
  vec
})
elements_up_df <- data.frame(padded_list_up)
write.csv(elements_up_df, file = "E:/paper-files/total_lcm_venn_elements_up_df.csv", 
          row.names = FALSE)

#repeat for downregulated genes
vectors_list_down <- list(sxcd_down, sbt_down, adeno_down, hgsoc_down)
elements_down <- calculate.overlap(vectors_list_down)
# Find the maximum length of the vectors
max_length_down <- max(sapply(elements_down, length))
# Pad the vectors with NA to make them of equal length
padded_list_down <- lapply(elements_down, function(vec) {
  length(vec) <- max_length_down
  vec
})
elements_down_df <- data.frame(padded_list_down)

write.csv(elements_down_df, file = "E:/paper-files/total_lcm_venn_elements_down_df.csv", row.names = FALSE)

################################import logcpm+c files##############################
logcpm_mouse <- read.csv('E:/paper-files/mlcm_total_logcpmc.csv', sep=',', header = TRUE)
logcpm_human <- read.csv('E:/paper-files/hlcm_total_logcpmc.csv', sep=',', header = TRUE)

#add gene symbols
logcpm_mouse <- as.data.frame(logcpm_mouse)
colnames(logcpm_mouse)[1] <- "ensgene"
logcpm_mouse <- rownames_to_column(logcpm_mouse)

logcpm_human <- as.data.frame(logcpm_human)
colnames(logcpm_human)[1] <- "ensgene"
logcpm_human <- rownames_to_column(logcpm_human)

logcpm_mouse <- left_join(x = logcpm_mouse,
                          y = grcm38 [, c("ensgene", "symbol", "entrez", "biotype", "description")], 
                          by = "ensgene")
logcpm_mouse <- logcpm_mouse %>%
  mutate(symbol = tolower(symbol))


logcpm_human <- left_join(x = logcpm_human,
                          y = grch38 [, c("ensgene", "symbol", "entrez", "biotype", "description")], 
                          by = "ensgene")
logcpm_human <- logcpm_human %>%
  mutate(symbol = tolower(symbol))

##################################filter logCPM data to common genes#########################
common_degs_up  <- na.omit(elements_up_df$a6) #filter for the a6 column which is the intersection
common_degs_up <- c(common_degs_up) #concatenate into vector

common_degs_down <- na.omit(elements_down_df$a6)
common_degs_down <- c(common_degs_down)

mouse_common_up <- logcpm_mouse %>%
  filter(symbol %in% common_degs_up)
mouse_common_down <- logcpm_mouse %>%
  filter(symbol %in% common_degs_down)

common_degs <- c(common_degs_up, common_degs_down)

mouse_common <- dplyr::bind_rows(mouse_common_up, mouse_common_down)

human_common_up <- logcpm_human %>%
  filter(symbol %in% common_degs_up)
human_common_down <- logcpm_human %>%
  filter(symbol %in% common_degs_down)

human_common <- dplyr::bind_rows(human_common_up, human_common_down)

mouse_names <- c("ensgene", "Adenoma_1", "Adenoma_2", "Adenoma_3", "Adenoma_4", "Adenoma_5", "Adenoma_6", "Adenoma_7", 
           "m_Control_1", "m_Control_2", "m_Control_3", "m_Control_4", "m_Control_5", "m_Control_6", "m_Control_7", 
           "Sex-cords_1", "Sex-cords_2", "Sex-cords_3", "Sex-cords_4", "Sex-cords_5", "Sex-cords_6", "Sex-cords_7", 
           "m_Control_8", "m_Control_9", "m_Control_10", "m_Control_11", "m_Control_12", "m_Control_13", "m_Control_14", 
           "m_Control_15", "m_Control_16", "m_Control_17", "m_Control_18", "m_Control_19", "m_Control_20", "m_Control_21", 
           "m_Control_22", "m_Control_23", "m_Control_24", "m_Control_25", "m_Control_26", "m_Control_27", "m_Control_28", 
           "symbol", "entrez", "biotype", "description")
human_names <- c("ensgene", "h_Control_1", "h_Control_2", "h_Control_3", "h_Control_4", "h_Control_5", 
                 "HGSOC_1", "HGSOC_2", "HGSOC_3", 
                 "h_Control_6", "h_Control_7", "h_Control_8", "h_Control_9", 
                 "SBT_1", "SBT_2", "SBT_3", "SBT_4", "SBT_5", "SBT_6", 
                 "symbol", "entrez", "biotype", "description")

mouse_common <- mouse_common[,-1]
human_common <- human_common[,-1]

colnames(mouse_common) <- mouse_names
colnames(human_common) <- human_names

mouse_common <- arrange(mouse_common, symbol)
human_common <- arrange(human_common, symbol)

common_symbol <- mouse_common$symbol

common_symbol <- toTitleCase(common_symbol)

#left with mouse_common and human_common


##there's a duplicate in the human dataframe, need to delete that##
# Check if the number of distinct symbols is equal to the total number of symbols
mouse_unique <- nrow(mouse_common) == nrow(distinct(mouse_common, symbol))
human_unique <- nrow(human_common) == nrow(distinct(human_common, symbol))

# Print the results
mouse_unique
human_unique

human_common <- human_common[-233,]

mouse_exp <- mouse_common[,2:43]
human_exp <- human_common[,2:19]

m_h_exp <- dplyr::bind_cols(mouse_exp, human_exp)

##add back the symbols as the first column

m_h_exp <- data.frame(Gene = common_symbol, m_h_exp)
rownames(m_h_exp) <- m_h_exp$Gene
m_h_exp <- m_h_exp[,2:61]

#left with a matrix that has gene symbols in alpha order for all of the samples
#the values of 
#####################################HEATMAP####################################
#make dataframe into matrix
m_h_exp <- as.matrix(m_h_exp)
m_h_exp_df <- as.data.frame(m_h_exp) #this only gives each individual log2CPM value

#so it looks like I filtered all the original files by the common_degs list (up and down)
#and then included the symbol, log2FC and padj values into a matrix.
###I have to figure out what all this means....I dont remember what the heck I was up to.

sxcd_lfc <- sxcd %>%
  filter(symbol %in% common_degs)%>%
  arrange(symbol)
sxcd_lfc <- as.matrix(dplyr::select(sxcd_lfc, symbol, log2FoldChange, padj))
adeno_lfc <- adeno %>%
  filter(symbol %in% common_degs)%>%
  arrange(symbol)
adeno_lfc <- dplyr::select(adeno_lfc, symbol, log2FoldChange, padj)
sbt_lfc <- sbt %>%
  filter(symbol %in% common_degs)%>%
  arrange(symbol)
sbt_lfc <- dplyr::select(sbt_lfc, symbol, log2FoldChange, padj)
hgsoc_lfc <- hgsoc %>%
  filter(symbol %in% common_degs) %>%
  arrange(symbol)
hgsoc_lfc <- dplyr::select(hgsoc_lfc, symbol, log2FoldChange, padj)

#tmprss2 is repeated so remove that gene

hgsoc_lfc <- hgsoc_lfc[-233,]
sbt_lfc <- sbt_lfc[-233,]

lfc_all <- dplyr::bind_cols(sxcd_lfc, adeno_lfc$log2FoldChange, adeno_lfc$padj, sbt_lfc$log2FoldChange, 
                            sbt_lfc$padj, hgsoc_lfc$log2FoldChange, hgsoc_lfc$padj)

new_col_names <- c("symbol", "Sexcords", "sc_padj", "Adenoma", "adeno_padj", "SBT", "sbt_padj", 
                   "HGSOC", "hg_padj")
colnames(lfc_all) <- new_col_names

sample_groups <- c("Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma", 
                   "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", 
                   "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", 
                   "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", 
                   "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", 
                   "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", 
                   "Human Control", "Human Control", "Human Control", "Human Control", "HGSOC", "HGSOC", "HGSOC", 
                   "Human Control", "Human Control", "Human Control", "Human Control", 
                   "SBT", "SBT", "SBT", "SBT", "SBT", "SBT")

sample_groups <- factor(sample_groups, levels = c("Mouse Control", "Sex cords", "Adenoma", 
                                                  "Human Control", "SBT", "HGSOC"))

anno_colours <- c("Mouse Control" = "gold", "Sex cords" = "plum3", "Adenoma" = "indianred", 
                  "Human Control" = "lightgreen", "SBT" = "cadetblue", "HGSOC" = "pink")

sample_groups <- factor(sample_groups, levels = names(anno_colours))

##make a dotplot of top 20 up and down genes
##need to make up and down regulated gene list...

top_dfs <- list(sxcd_lfc = sxcd_lfc, 
                adeno_lfc = adeno_lfc, 
                sbt_lfc = sbt_lfc, 
                hgsoc_lfc = hgsoc_lfc)

# Loop over each dataframe in the list
for (name in names(top_dfs)) {
  df <- as.data.frame(top_dfs[[name]])
  df$log2FoldChange <- as.numeric(df$log2FoldChange)
  # Sort by log2FoldChange in decreasing order for upregulated genes
  df_sorted_up <- df[order(df$log2FoldChange, decreasing = TRUE), ]
  top20_up <- head(df_sorted_up, 20)
    # Sort by log2FoldChange in increasing order for downregulated genes
  df_sorted_down <- df[order(df$log2FoldChange, decreasing = FALSE), ]
  top20_down <- head(df_sorted_down, 20)
    # Store the top 20 dataframes in the global environment
  assign(paste0(name, "_top20_up"), top20_up)
  assign(paste0(name, "_top20_down"), top20_down)
}

# Access the results they now should be of the structure sxcd_lfc_top20_up
#add the phenotype to each dataframe
sxcd_lfc_top20_up$condition <- "Sexcords"
adeno_lfc_top20_up$condition <- "Adenoma"
sbt_lfc_top20_up$condition <- "SBT"
hgsoc_lfc_top20_up$condition <- "HGSOC"
sxcd_lfc_top20_down$condition <- "Sexcords"
adeno_lfc_top20_down$condition <- "Adenoma"
sbt_lfc_top20_down$condition <- "SBT"
hgsoc_lfc_top20_down$condition <- "HGSOC"

sxcd_lfc_top20_up$regulation <- "Up"
adeno_lfc_top20_up$regulation <- "Up"
sbt_lfc_top20_up$regulation <- "Up"
hgsoc_lfc_top20_up$regulation <- "Up"
sxcd_lfc_top20_down$regulation <- "Down"
adeno_lfc_top20_down$regulation <- "Down"
sbt_lfc_top20_down$regulation <- "Down"
hgsoc_lfc_top20_down$regulation <- "Down"

# List of dataframes
dataframes_list <- list(sxcd_lfc_top20_up, 
                        adeno_lfc_top20_up, 
                        sbt_lfc_top20_up, 
                        hgsoc_lfc_top20_up, 
                        sxcd_lfc_top20_down, 
                        adeno_lfc_top20_down, 
                        sbt_lfc_top20_down, 
                        hgsoc_lfc_top20_down)

top_genes_all <- do.call(rbind, dataframes_list)

levels_condition <- c("Sexcords", "Adenoma", "SBT", "HGSOC")
levels_regulation <- c("Up", "Down")

# Relevel condition column
top_genes_all$condition <- factor(top_genes_all$condition, levels = levels_condition)

# Relevel regulation column
top_genes_all$regulation <- factor(top_genes_all$regulation, levels = levels_regulation)

# Print the updated dataframe
print(top_genes_all)

#convert padj to numeric to calculate -log for the dotplot
top_genes_all$padj <- as.numeric(as.character(top_genes_all$padj))
top_genes_all$log_padj <- -log10(top_genes_all$padj)
top_genes_all$symbol <- toTitleCase(top_genes_all$symbol)

top_genes_up <- top_genes_all[top_genes_all$regulation == "Up", ]
top_genes_down <- top_genes_all[top_genes_all$regulation == "Down", ]

up <- ggplot(top_genes_up, aes(x = condition, y = symbol, size = log2FoldChange, color = log_padj)) +
  geom_point() +  # Separate points by regulation (Up/Down)
  scale_color_gradient(low = "blue", high = "red") +  # Color gradient for padj
  scale_size_continuous(range = c(1, 10)) +  # Adjust size range as needed
  labs(x = "Phenotype", y = "Genes", size = "log2FC", color = "-log10padj") +  
  theme_grey()  # Optional: Adjust theme as needed

down <- ggplot(top_genes_down, aes(x = condition, y = symbol, size = log2FoldChange, color = log_padj)) +
  geom_point() +  # Separate points by regulation (Up/Down)
  scale_color_gradient(low = "blue", high = "red") +  # Color gradient for padj
  scale_size_continuous(range = c(1, 10)) +  # Adjust size range as needed
  labs(x = "Phenotype", y = "Genes", size = "log2FC", color = "-log10padj") +  
  theme_grey()  # Optional: Adjust theme as needed

# Display the plot
print(up)
print(down)

## not sure I like that representation of the data

##maybe too many genes to represent with gene expression data...
##maybe correlation matrix instead?

corr_value <- cor(m_h_exp)
corrplot <- pheatmap(corr_value, 
                     fontsize_row = 8, 
                     fontsize_col = 8, 
                     main = "Correlation Values of Sample Transcriptomes")
print(corrplot)

###get prepped for complex heatmap
group <- c("Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma", 
                   "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", 
                   "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", 
                   "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", 
                   "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", 
                   "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", "Mouse Control", 
                   "Human Control", "Human Control", "Human Control", "Human Control", "Human Control", "HGSOC", "HGSOC", "HGSOC", 
                   "Human Control", "Human Control", "Human Control", "Human Control", 
                   "SBT", "SBT", "SBT", "SBT", "SBT", "SBT")
condition <- c("Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma", "Adenoma", 
                       "KO Stroma", "KO Stroma", "KO Stroma", "KO Stroma", "KO Stroma", "KO Stroma", "KO Stroma", 
                       "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", "Sex cords", 
                       "Mature GCs", "Mature GCs", "Mature GCs", "Mature GCs", "Mature GCs", "Mature GCs", "Mature GCs", 
                       "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs", "Primary GCs", 
                       "WT Stroma", "WT Stroma", "WT Stroma", "WT Stroma", "WT Stroma", "WT Stroma", "WT Stroma", 
                       "Benign", "Benign", "Benign", "Benign", "Benign", "HGSOC", "HGSOC", "HGSOC", 
                       "Normal", "Normal", "Normal", "Normal", 
                       "SBT", "SBT", "SBT", "SBT", "SBT", "SBT")
species <- c("Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", 
                    "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", 
                    "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", 
                    "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", 
                    "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", 
                    "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", "Mouse", 
                    "Human",  "Human",  "Human",  "Human",  "Human",  "Human",  "Human",  "Human", 
                    "Human",  "Human",  "Human",  "Human", 
                    "Human",  "Human",  "Human",  "Human",  "Human",  "Human")

metadata <- data.frame(species, condition, group)
rownames(metadata) <- colnames(m_h_exp)
all(rownames(metadata) == colnames(m_h_exp))

order_condition <- c("Adenoma", "Sex cords", "Primary GCs", "Mature GCs", "WT Stroma", "KO Stroma", 
                     "HGSOC", "SBT", "Normal", "Benign")
order_group <- c("Mouse Control", "Sex cords", "Adenoma", 
                 "Human Control", "SBT", "HGSOC")
order_species <- c("Mouse", "Human")
                                  
#relevel the metadata to match the PCA plot
metadata$condition <- factor(metadata$condition, levels = order_condition)
metadata$group <- factor(metadata$group, levels = order_group)
metadata$species <- factor(metadata$species, levels = order_species)

condition_colours <- c("Adenoma" = "#D81B60", "Mature GCs" = "#1E88E5", 
                       "Primary GCs" = "#5D286B", "Sex cords" = "#004D40", 
                       "WT Stroma" = "#FFC107", "KO Stroma" = "darkorange", 
                       "HGSOC" = "#D81B60", "Normal" = "#1E88E5", "SBT" = "#004D40", "Benign" = "#5D286B")
group_colours <- c("Mouse Control" = "#FFC107", "Sex cords" = "#004D40", "Adenoma" = "#D81B60", 
                   "Human Control" = "darkorange", "SBT" = "magenta", "HGSOC" = "seagreen")
species_colours <- c("Mouse" = "black", "Human" = "honeydew4")

#colour blind safe

ha_top <- HeatmapAnnotation(
  species = metadata$species, 
  condition = metadata$condition,
  group = metadata$group,
  col = list(species = species_colours, 
             condition = condition_colours, 
            group = group_colours),
  annotation_height = unit(3, "mm"), 
  annotation_width = unit(0.5, 'cm'), 
  gap = unit(0.5, 'mm'), 
  border = TRUE, 
  annotation_legend_param = list(
    species = list(
      nrow = 2, 
      title = "Species", 
      title_position = 'topleft',
      legend_direction = 'horizontal',
      title_gp = gpar(fontsize =12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'plain')),
    condition = list(
      nrow = 10, 
      title = "Condition", 
      title_position = 'topleft',
      legend_direction = 'horizontal',
      title_gp = gpar(fontsize =12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'plain')),
    group = list(
      nrow = 6,
      title = 'Group',
      title_position = 'topleft',
      legend_direction = 'horizontal',
      title_gp = gpar(fontsize = 12, fontface = 'bold'),
      labels_gp = gpar(fontsize = 12, fontface = 'plain'))))

#row annotation (order of variables is th eorder of the annotation bars)
ha_row <- HeatmapAnnotation(which = "row",
                            species = metadata$species,                          
                            condition = metadata$condition,
                            group = metadata$group, 
                            col = list(species = species_colours, 
                                       condition = condition_colours, 
                                       group = group_colours
                            ),
                            annotation_height = 0.3, 
                            annotation_width = unit(1, 'cm'), 
                            gap = unit(1, 'mm'), 
                            border = TRUE, 
                            show_legend = FALSE, 
                            show_annotation_name = FALSE, 
                            annotation_legend_param = list(
                              species = list(
                                nrow = 2,
                                title = 'Species',
                                title_position = 'topcenter',
                                legend_direction = 'vertical',
                                title_gp = gpar(fontsize = 12, fontface = 'plain'),
                                labels_gp = gpar(fontsize = 12, fontface = 'plain')), 
                              condition = list(
                                nrow = 10, 
                                title = "Condition", 
                                title_position = 'topleft',
                                legend_direction = 'horizontal',
                                title_gp = gpar(fontsize =12, fontface = 'bold'),
                                labels_gp = gpar(fontsize = 12, fontface = 'plain')),
                              group = list(
                                nrow = 6,
                                title = 'Group',
                                title_position = 'topleft',
                                legend_direction = 'horizontal',
                                title_gp = gpar(fontsize = 12, fontface = 'bold'),
                                labels_gp = gpar(fontsize = 12, fontface = 'plain'))))

corrplot2 <- ComplexHeatmap::Heatmap(corr_value, 
                                     name = "Correlation of Mouse and Human Common DEGs", 
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






sxcd_lfc <- as.data.frame(sxcd_lfc)

#try heat map with logFC >= 2
logfc_dfs <- list(sxcd_lfc = sxcd_lfc, 
                  adeno_lfc = adeno_lfc, 
                  sbt_lfc = sbt_lfc, 
                  hgsoc_lfc =  hgsoc_lfc)

for (name in names(logfc_dfs)) {
  df <- logfc_dfs[[name]]
  df$log2FoldChange <- as.numeric(df$log2FoldChange)
  # Filter out NA values and apply the filter condition
  filtered_df <- df %>%
    filter(log2FoldChange >= 2)
    # Assign the filtered dataframe to a new variable dynamically
  assign(paste0(name, "_log2"), filtered_df)
}

# Print the filtered dataframes to check
print(sxcd_lfc_log2)
print(adeno_lfc_log2)
print(sbt_lfc_log2)
print(hgsoc_lfc_log2)

sxcd_lfc$log2FoldChange <- as.numeric(sxcd_lfc$log2FoldChange)
sxcd_lfc_2 <- sxcd_lfc %>%
  filter(log2FoldChange >= 2)


sxcd_lfc_2 <- as.matrix(filter(sxcd_lfc, log2FoldChange >= 2))
sxcd_lfc <- as.matrix(dplyr::select(sxcd_lfc, symbol, log2FoldChange, padj))
adeno_lfc <- adeno %>%
  filter(symbol %in% common_degs)%>%
  arrange(symbol)
adeno_lfc <- dplyr::select(adeno_lfc, symbol, log2FoldChange, padj)
sbt_lfc <- sbt %>%
  filter(symbol %in% common_degs)%>%
  arrange(symbol)
sbt_lfc <- dplyr::select(sbt_lfc, symbol, log2FoldChange, padj)
hgsoc_lfc <- hgsoc %>%
  filter(symbol %in% common_degs) %>%
  arrange(symbol)
hgsoc_lfc <- dplyr::select(hgsoc_lfc, symbol, log2FoldChange, padj)
heatmap <- Heatmap(lfc_all$Sexcords, cluster_rows = TRUE, cluster_columns = TRUE, 
                              name = "log2CPM + c", 
                              rect_gp = gpar(col = "white", lwd = 0.25), 
                              clustering_distance_rows = "euclidean", 
                              row_names_gp = gpar(fontsize = 14), 
                              heatmap_legend_param = list(title_gp = gpar(fontsize = 14), 
                                                          labels_gp = gpar(fontsize = 14)))   
print(heatmap)


#heat_m_h_common <- columnAnnotation(heatmap_m_h_common, 
                                 #   col = list(SampleGroup = sample_groups), 
                                  #  annotation_legend_param = list(SampleGroup = list(title = "Groups")))
col_anno <- columnAnnotation(Group = sample_groups, 
                             col = list(Group = anno_colours), 
                             border = TRUE)
heatmap_lfc <- Heatmap(sxcd_lfc)
heatmap_m_h_common <- heatmap_m_h_common %v% col_anno
heatmap_m_h_common <- heatmap_m_h_common + heatmap_lfc
                             
draw(heatmap_m_h_common)
draw(heatmap_lfc)


#################################

###############################ADENO GO#########################################
##make background genes list with correct gene symbol case.
logcpm_mouse_go <- read.csv('E:/paper-files/mlcm_total_logcpmc.csv', sep=',', header = TRUE)
logcpm_human_go <- read.csv('E:/paper-files/hlcm_total_logcpmc.csv', sep=',', header = TRUE)

#add gene symbols
logcpm_mouse_go <- as.data.frame(logcpm_mouse_go)
colnames(logcpm_mouse_go)[1] <- "ensgene"

logcpm_human_go <- as.data.frame(logcpm_human_go)
colnames(logcpm_human_go)[1] <- "ensgene"

logcpm_mouse_go <- left_join(x = logcpm_mouse_go,
                          y = grcm38 [, c("ensgene", "symbol", "entrez", "biotype", "description")], 
                          by = "ensgene")

logcpm_human_go <- left_join(x = logcpm_human_go,
                          y = grch38 [, c("ensgene", "symbol", "entrez", "biotype", "description")], 
                          by = "ensgene")

mouse_background <- logcpm_mouse_go %>%
  filter(biotype == "protein_coding") %>%
  dplyr::select(symbol)

human_background <- logcpm_human_go %>%
  filter(biotype == "protein_coding") %>%
  dplyr::select(symbol)

human_background <- unique(human_background) #get rid of gene duplicates
human_background <- human_background[human_background != ""] #get rid of blanks

#get gene lists again but in correct case
sxcd_go <- read.csv('E:/paper-files/mlcm_total_sxcd_v_control_deseq_xlfcshrink.csv', sep=',', header = TRUE)
adeno_go <- read.csv('E:/paper-files/mlcm_total_adeno_v_control_deseq_xlfcshrink.csv', sep=',', header = TRUE)
sbt_go <- read.csv('E:/paper-files/hlcm_total_SBT_v_control_deseq_xlfcshrink.csv', sep=',', header = TRUE)
hgsoc_go <- read.csv('E:/paper-files/hlcm_total_HGSOC_v_control_deseq_xlfcshrink.csv', sep=',', header = TRUE)

#filter each dataset for mRNAs
sxcd_go_up <- filter(sxcd_go, log2FoldChange >= 0, biotype == "protein_coding", symbol != "")
sxcd_go_down <- filter(sxcd_go, log2FoldChange < 0, biotype == "protein_coding", symbol != "")
adeno_go_up <- filter(adeno_go, log2FoldChange >= 0, biotype == "protein_coding", symbol != "")
adeno_go_down <- filter(adeno_go, log2FoldChange < 0, biotype == "protein_coding", symbol != "")
sbt_go_up <- filter(sbt_go, log2FoldChange >= 0, biotype == "protein_coding", symbol != "")
sbt_go_down <- filter(sbt_go, log2FoldChange < 0, biotype == "protein_coding", symbol != "")
hgsoc_go_up <- filter(hgsoc_go, log2FoldChange >= 0, biotype == "protein_coding", symbol != "")
hgsoc_go_down <- filter(hgsoc_go, log2FoldChange < 0, biotype == "protein_coding", symbol != "")

sxcd_go_up <- sxcd_go_up$symbol
sxcd_go_down <- sxcd_go_down$symbol
adeno_go_up <- adeno_go_up$symbol
adeno_go_down <- adeno_go_down$symbol
sbt_go_up <- sbt_go_up$symbol
sbt_go_down <- sbt_go_down$symbol
hgsoc_go_up <- hgsoc_go_up$symbol
hgsoc_go_down <- hgsoc_go_down$symbol

#test <- write.csv(sbt_go_up, 'E:/paper-files/test_sbt_up.csv', row.names = FALSE)

GO_adeno_up <- enrichGO(gene = adeno_go_up, OrgDb = "org.Mm.eg.db", 
                       keyType = "SYMBOL", ont = "ALL", 
                       pAdjustMethod = "fdr", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05, 
                       minGSSize = 15, 
                       maxGSSize = 500, 
                       universe = mouse_background)

GO_adeno_up <- as.data.frame(GO_adeno_up)
GO_adeno_up['direction'] = "Up"

GO_adeno_down <- enrichGO(gene = adeno_go_down, OrgDb = "org.Mm.eg.db", 
                         keyType = "SYMBOL", ont = "ALL", 
                         pAdjustMethod = "fdr", 
                         pvalueCutoff = 0.05, 
                         qvalueCutoff = 0.05, 
                         minGSSize = 15, 
                         maxGSSize = 500, 
                         universe = mouse_background)

GO_adeno_down <- as.data.frame(GO_adeno_down)
GO_adeno_down['direction'] = "Down"

GO_adeno <- dplyr::bind_rows(GO_adeno_up, GO_adeno_down)
write.csv(GO_adeno, "E:/paper-files/GO_adeno.csv", row.names = TRUE)
#######################################SXCD GO##################################
GO_sxcd_up <- enrichGO(gene = sxcd_go_up, OrgDb = "org.Mm.eg.db", 
                        keyType = "SYMBOL", ont = "ALL", 
                        pAdjustMethod = "fdr", 
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.05, 
                        minGSSize = 15, 
                        maxGSSize = 500, 
                        universe = mouse_background)

GO_sxcd_up <- as.data.frame(GO_sxcd_up)
GO_sxcd_up['direction'] = "Up"

GO_sxcd_down <- enrichGO(gene = sxcd_go_down, OrgDb = "org.Mm.eg.db", 
                          keyType = "SYMBOL", ont = "ALL", 
                          pAdjustMethod = "fdr", 
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.05, 
                          minGSSize = 15, 
                          maxGSSize = 500, 
                          universe = mouse_background)

GO_sxcd_down <- as.data.frame(GO_sxcd_down)
GO_sxcd_down['direction'] = "Down"

GO_sxcd <- dplyr::bind_rows(GO_sxcd_up, GO_sxcd_down)
write.csv(GO_sxcd, "E:/paper-files/GO_sxcd.csv", row.names = TRUE)

########################################SBT GO##################################
GO_sbt_up <- enrichGO(gene = sbt_go_up, OrgDb = "org.Hs.eg.db", 
                            keyType = "SYMBOL", ont = "ALL", 
                      pAdjustMethod = "fdr", 
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.05, 
                      minGSSize = 15, 
                      maxGSSize = 500, 
                      universe = human_background)

GO_sbt_up <- as.data.frame(GO_sbt_up)
GO_sbt_up['direction'] = "Up"

#sbt_go_down <- sbt_go_down[sbt_go_down != "" & !is.na(sbt_go_down)]
#sbt_go_down <- unique(sbt_go_down)

#for some reason the sbt_down doesn't work when ontology is set to all
#do individually and then merge
#GO_sbt_down <- enrichGO(gene = sbt_go_down, OrgDb = "org.Hs.eg.db", 
 #                     keyType = "SYMBOL", ont = "ALL", 
  #                    pAdjustMethod = "fdr", 
   #                   pvalueCutoff = 0.05, 
    #                  qvalueCutoff = 0.05, 
     #                 minGSSize = 15, 
      #                maxGSSize = 500, 
       #               universe = human_background)
GO_sbt_down_bp <- enrichGO(gene = sbt_go_down, OrgDb = "org.Hs.eg.db", 
                           keyType = "SYMBOL", ont = "BP", 
                           pAdjustMethod = "fdr", 
                           pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.05, 
                           minGSSize = 15, 
                           maxGSSize = 500, 
                           universe = human_background)

GO_sbt_down_bp <- as.data.frame(GO_sbt_down_bp)
GO_sbt_down_bp <- GO_sbt_down_bp %>%
  mutate(ONTOLOGY = "BP") %>%
  dplyr::select(ONTOLOGY, everything())

GO_sbt_down_mf <- enrichGO(gene = sbt_go_down, OrgDb = "org.Hs.eg.db", 
                           keyType = "SYMBOL", ont = "MF", 
                           pAdjustMethod = "fdr", 
                           pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.05, 
                           minGSSize = 15, 
                           maxGSSize = 500, 
                           universe = human_background)

GO_sbt_down_mf <- as.data.frame(GO_sbt_down_mf)
GO_sbt_down_mf <- GO_sbt_down_mf %>%
  mutate(ONTOLOGY = "MF") %>%
  dplyr::select(ONTOLOGY, everything())

GO_sbt_down_cc <- enrichGO(gene = sbt_go_down, OrgDb = "org.Hs.eg.db", 
                           keyType = "SYMBOL", ont = "CC", 
                           pAdjustMethod = "fdr", 
                           pvalueCutoff = 0.05, 
                           qvalueCutoff = 0.05, 
                           minGSSize = 15, 
                           maxGSSize = 500, 
                           universe = human_background)

GO_sbt_down_cc <- as.data.frame(GO_sbt_down_cc)
GO_sbt_down_cc <- GO_sbt_down_cc %>%
  mutate(ONTOLOGY = "CC") %>%
  dplyr::select(ONTOLOGY, everything())

GO_sbt_down_bp <- as.data.frame(GO_sbt_down_bp)
GO_sbt_down_bp['direction'] = "Down"
GO_sbt_down_mf <- as.data.frame(GO_sbt_down_mf)
GO_sbt_down_mf['direction'] = "Down"
GO_sbt_down_cc <- as.data.frame(GO_sbt_down_cc)
GO_sbt_down_cc['direction'] = "Down"

GO_sbt <- dplyr::bind_rows(GO_sbt_up, GO_sbt_down_bp, GO_sbt_down_mf, GO_sbt_down_cc)
write.csv(GO_sbt, "E:/paper-files/GO_sbt.csv", row.names = TRUE)

########################################HGSOC GO################################
GO_hgsoc_up <- enrichGO(gene = hgsoc_go_up, OrgDb = "org.Hs.eg.db", 
                       keyType = "SYMBOL", ont = "ALL", 
                       pAdjustMethod = "fdr", 
                       pvalueCutoff = 0.05, 
                       qvalueCutoff = 0.05, 
                       minGSSize = 15, 
                       maxGSSize = 500, 
                       universe = human_background)

GO_hgsoc_up <- as.data.frame(GO_hgsoc_up)
GO_hgsoc_up['direction'] = "Up"

GO_hgsoc_down <- enrichGO(gene = hgsoc_go_down, OrgDb = "org.Hs.eg.db", 
                        keyType = "SYMBOL", ont = "ALL", 
                        pAdjustMethod = "fdr", 
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.05, 
                        minGSSize = 15, 
                        maxGSSize = 500, 
                        universe = human_background)

GO_hgsoc_down <- as.data.frame(GO_hgsoc_down)
GO_hgsoc_down['direction'] = "Down"

GO_hgsoc <- dplyr::bind_rows(GO_hgsoc_up, GO_hgsoc_down)
write.csv(GO_hgsoc, "E:/paper-files/GO_hgsoc.csv", row.names = TRUE)

#combine all GO mouse and human into one data frame and then export

GO_sxcd['group'] = "Sex cords"
GO_sxcd['species'] = "mouse"
GO_adeno['group'] = "Adenoma"
GO_adeno['species'] = "mouse"
GO_sbt['group'] = "SBT"
GO_sbt['species'] = "human"
GO_hgsoc['group'] = "HGSOC"
GO_hgsoc['species'] = "human"

GO_mrna <- dplyr::bind_rows(GO_sxcd, GO_adeno, GO_sbt, GO_hgsoc)
write.csv(GO_mrna, "E:/paper-files/GO_mrna.csv", row.names = TRUE)
##################Venn intersection for each GO and Regulation Direction########
#####################################UPREGULATED GO VENN########################
#Adenoma UP
GO_adeno_up_venn_bp <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "BP", direction == "Up", group == "Adenoma")
GO_adeno_up_venn_bp <- GO_adeno_up_venn_bp$ID

GO_adeno_up_venn_cc <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "CC", direction == "Up", group == "Adenoma")
GO_adeno_up_venn_cc <- GO_adeno_up_venn_cc$ID

GO_adeno_up_venn_mf <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "MF", direction == "Up", group == "Adenoma")
GO_adeno_up_venn_mf <- GO_adeno_up_venn_mf$ID

#Sxcd Up
GO_sxcd_up_venn_bp <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "BP" & direction == "Up", group == "Sex cords")
GO_sxcd_up_venn_bp <- GO_sxcd_up_venn_bp$ID

GO_sxcd_up_venn_cc <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "CC" & direction == "Up", group == "Sex cords")
GO_sxcd_up_venn_cc <- GO_sxcd_up_venn_cc$ID

GO_sxcd_up_venn_mf <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "MF" & direction == "Up", group == "Sex cords")
GO_sxcd_up_venn_mf <- GO_sxcd_up_venn_mf$ID
#HGSOC Up
GO_hgsoc_up_venn_bp <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "BP" & direction == "Up", group == "HGSOC")
GO_hgsoc_up_venn_bp <- GO_hgsoc_up_venn_bp$ID

GO_hgsoc_up_venn_cc <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "CC" & direction == "Up", group == "HGSOC")
GO_hgsoc_up_venn_cc <- GO_hgsoc_up_venn_cc$ID

GO_hgsoc_up_venn_mf <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "MF" & direction == "Up", group == "HGSOC")
GO_hgsoc_up_venn_mf <- GO_hgsoc_up_venn_mf$ID

#SBT
GO_sbt_up_venn_bp <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "BP" & direction == "Up", group == "SBT")
GO_sbt_up_venn_bp <- GO_sbt_up_venn_bp$ID

GO_sbt_up_venn_cc <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "CC" & direction == "Up", group == "SBT")
GO_sbt_up_venn_cc <- GO_sbt_up_venn_cc$ID

GO_sbt_up_venn_mf <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "MF" & direction == "Up", group == "SBT")
GO_sbt_up_venn_mf <- GO_sbt_up_venn_mf$ID

up_bp <- list(Sexcords = GO_sxcd_up_venn_bp, 
              SBT = GO_sbt_up_venn_bp, 
              Adenoma = GO_adeno_up_venn_bp, 
              HGSOC = GO_hgsoc_up_venn_bp)

venn_up_bp <- venn.diagram(x = up_bp,
                           category.names = c("Sexcords", "SBT", "Adenoma", "HGSOC"), 
                           filename = "E:/paper-files/images/go_bp_venn_up_mh_id.tiff",
                           disable.logging = FALSE, 
                           imagetype = "tiff", 
                           main = "Mouse versus Human - Biological Processes", 
                           main.fontface = "bold", 
                           sub = "Upregulated DEGs", 
                           col = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                           fill = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                           alpha = 0.3, 
                           fontfamily = "sans", 
                           resolution = 800, 
                           width = 6, 
                           height = 6, 
                           units = "in", 
                           sub.fontface = "bold", 
                           cat.fontfamily = "sans", 
                           cat.fontface = "bold", 
                           cex = 3,
                           cat.cex = 1.5)

up_cc <- list(Sexcords = GO_sxcd_up_venn_cc, 
              SBT = GO_sbt_up_venn_cc, 
              Adenoma = GO_adeno_up_venn_cc, 
              HGSOC = GO_hgsoc_up_venn_cc)

venn_up_cc <- venn.diagram(x = up_cc,
                           category.names = c("Sexcords", "SBT", "Adenoma", "HGSOC"), 
                           filename = "E:/paper-files/images/go_cc_venn_up_mh_id.tiff",
                           disable.logging = FALSE, 
                           imagetype = "tiff", 
                           main = "Mouse versus Human - Cellular Components", 
                           main.fontface = "bold", 
                           sub = "Upregulated DEGs", 
                           col = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                           fill = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                           alpha = 0.3, 
                           fontfamily = "sans", 
                           resolution = 800, 
                           width = 6, 
                           height = 6, 
                           units = "in", 
                           sub.fontface = "bold", 
                           cat.fontfamily = "sans", 
                           cat.fontface = "bold", 
                           cex = 3,
                           cat.cex = 1.5)

up_mf <- list(Sexcords = GO_sxcd_up_venn_mf, 
              SBT = GO_sbt_up_venn_mf, 
              Adenoma = GO_adeno_up_venn_mf, 
              HGSOC = GO_hgsoc_up_venn_mf)

venn_up_mf <- venn.diagram(x = up_mf,
                           category.names = c("Sexcords", "SBT", "Adenoma", "HGSOC"), 
                           filename = "E:/paper-files/images/go_mf_venn_up_mh_id.tiff",
                           disable.logging = FALSE, 
                           imagetype = "tiff", 
                           main = "Mouse versus Human - Molecular Functions", 
                           main.fontface = "bold", 
                           sub = "Upregulated DEGs", 
                           col = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                           fill = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                           alpha = 0.3, 
                           fontfamily = "sans", 
                           resolution = 800, 
                           width = 6, 
                           height = 6, 
                           units = "in", 
                           sub.fontface = "bold", 
                           cat.fontfamily = "sans", 
                           cat.fontface = "bold", 
                           cex = 3,
                           cat.cex = 1.5)
###upset plots of common genes:
##sets created
#up_bp, down_bp, up_cc, down_cc, up_mf, down_mf

CairoPNG("E:/paper-files/images/up_bp_upset.png", width = 6, height = 6, units = "in", res = 800, 
         bg = "transparent")

upset(
  fromList(up_bp),
  order.by = "freq",  # Order by frequency
  point.size = 3.5,  # Size of points in the matrix
  line.size = 1,  # Size of lines connecting sets
  main.bar.color = "seagreen",  # Color of the main bar
  sets.bar.color = "gray",  # Color of the sets bar
  matrix.color = "red3",  # Color of the matrix points
  text.scale = c(1.8, 1.8, 1.2, 1.2, 1.8, 1.8),  # Scale of text elements
  sets.x.label = "No. of Terms in Set",  # Label for sets bar
  keep.order = TRUE,  # Keep the order of sets
  empty.intersections = "on"  # Show empty intersections
)
# Close the graphics device
dev.off()

CairoPNG("E:/paper-files/images/up_cc_upset.png", width = 6, height = 6, units = "in", res = 800, 
         bg = "transparent")
upset(
  fromList(up_cc),
  order.by = "freq",  # Order by frequency
  point.size = 3.5,  # Size of points in the matrix
  line.size = 1,  # Size of lines connecting sets
  main.bar.color = "darkorange2",  # Color of the main bar
  sets.bar.color = "gray",  # Color of the sets bar
  matrix.color = "red3",  # Color of the matrix points
  text.scale = c(1.8, 1.8, 1.4, 1.4, 1.8, 1.8),  # Scale of text elements
  sets.x.label = "No. of Terms in Set",  # Label for sets bar
  keep.order = TRUE,  # Keep the order of sets
  empty.intersections = "on"  # Show empty intersections
)
# Close the graphics device
dev.off()

CairoPNG("E:/paper-files/images/up_mf_upset.png", width = 6, height = 6, units = "in", res = 800, 
         bg = "transparent")
upset(
  fromList(up_mf),
  order.by = "freq",  # Order by frequency
  point.size = 3.5,  # Size of points in the matrix
  line.size = 1,  # Size of lines connecting sets
  main.bar.color = "mediumpurple",  # Color of the main bar
  sets.bar.color = "gray",  # Color of the sets bar
  matrix.color = "red3",  # Color of the matrix points
  text.scale = c(1.8, 1.8, 1.4, 1.4, 1.8, 1.8),  # Scale of text elements
  sets.x.label = "No. of Terms in Set",  # Label for sets bar
  keep.order = TRUE,  # Keep the order of sets
  empty.intersections = "on"  # Show empty intersections
)
# Close the graphics device
dev.off()
##############################DOWNREGULATED GO VENN#############################
#Adenoma down
GO_adeno_down_venn_bp <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "BP", direction == "Down", group == "Adenoma")
GO_adeno_down_venn_bp <- GO_adeno_down_venn_bp$ID

GO_adeno_down_venn_cc <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "CC", direction == "Down", group == "Adenoma")
GO_adeno_down_venn_cc <- GO_adeno_down_venn_cc$ID

GO_adeno_down_venn_mf <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "MF", direction == "Down", group == "Adenoma")
GO_adeno_down_venn_mf <- GO_adeno_down_venn_mf$ID

#Sxcd down
GO_sxcd_down_venn_bp <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "BP" & direction == "Down", group == "Sex cords")
GO_sxcd_down_venn_bp <- GO_sxcd_down_venn_bp$ID

GO_sxcd_down_venn_cc <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "CC" & direction == "Down", group == " Sex cords")
GO_sxcd_down_venn_cc <- GO_sxcd_down_venn_cc$ID

GO_sxcd_down_venn_mf <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "MF" & direction == "Down", group == "Sex cords")
GO_sxcd_down_venn_mf <- GO_sxcd_down_venn_mf$ID
#hgsoc down
GO_hgsoc_down_venn_bp <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "BP" & direction == "Down", group == "HGSOC")
GO_hgsoc_down_venn_bp <- GO_hgsoc_down_venn_bp$ID

GO_hgsoc_down_venn_cc <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "CC" & direction == "Down", group == "HGSOC")
GO_hgsoc_down_venn_cc <- GO_hgsoc_down_venn_cc$ID

GO_hgsoc_down_venn_mf <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "MF" & direction == "Down", group == "HGSOC")
GO_hgsoc_down_venn_mf <- GO_hgsoc_down_venn_mf$ID

#SBT down
GO_sbt_down_venn_bp <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "BP" & direction == "Down", group == "SBT")
GO_sbt_down_venn_bp <- GO_sbt_down_venn_bp$ID

GO_sbt_down_venn_cc <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "CC" & direction == "Down", group == "SBT")
GO_sbt_down_venn_cc <- GO_sbt_down_venn_cc$ID

GO_sbt_down_venn_mf <- GO_mrna %>%
  dplyr::select(ONTOLOGY, ID, direction, group) %>%
  filter(ONTOLOGY == "MF" & direction == "Down", group == "SBT")
GO_sbt_down_venn_mf <- GO_sbt_down_venn_mf$ID

down_bp <- list(Sexcords = GO_sxcd_down_venn_bp, 
              SBT = GO_sbt_down_venn_bp, 
              Adenoma = GO_adeno_down_venn_bp, 
              HGSOC = GO_hgsoc_down_venn_bp)

venn_up_bp <- venn.diagram(x = down_bp,
                           category.names = c("Sexcords", "SBT", "Adenoma", "HGSOC"), 
                           filename = "E:/paper-files/images/go_bp_venn_down_mh_id.tiff",
                           disable.logging = FALSE, 
                           imagetype = "tiff", 
                           main = "Mouse versus Human - Biological Processes", 
                           main.fontface = "bold", 
                           sub = "Downregulated GO", 
                           col = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                           fill = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                           alpha = 0.3, 
                           fontfamily = "sans", 
                           resolution = 800, 
                           width = 6, 
                           height = 6, 
                           units = "in", 
                           sub.fontface = "bold", 
                           cat.fontfamily = "sans", 
                           cat.fontface = "bold", 
                           cex = 3,
                           cat.cex = 1.5)

down_cc <- list(Sexcords = GO_sxcd_down_venn_cc, 
                SBT = GO_sbt_down_venn_cc, 
                Adenoma = GO_adeno_down_venn_cc, 
                HGSOC = GO_hgsoc_down_venn_cc)

venn_down_cc <- venn.diagram(x = down_cc,
                           category.names = c("Sexcords", "SBT", "Adenoma", "HGSOC"), 
                           filename = "E:/paper-files/images/go_cc_venn_down_mh_id.tiff",
                           disable.logging = FALSE, 
                           imagetype = "tiff", 
                           main = "Mouse versus Human - Cellular Components", 
                           main.fontface = "bold", 
                           sub = "Downregulated GO", 
                           col = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                           fill = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                           alpha = 0.3, 
                           fontfamily = "sans", 
                           resolution = 800, 
                           width = 6, 
                           height = 6, 
                           units = "in", 
                           sub.fontface = "bold", 
                           cat.fontfamily = "sans", 
                           cat.fontface = "bold", 
                           cex = 3,
                           cat.cex = 1.5)

down_mf <- list(Sexcords = GO_sxcd_down_venn_mf, 
                SBT = GO_sbt_down_venn_mf, 
                Adenoma = GO_adeno_down_venn_mf, 
                HGSOC = GO_hgsoc_down_venn_mf)

venn_down_mf <- venn.diagram(x = down_mf,
                             category.names = c("Sexcords", "SBT", "Adenoma", "HGSOC"), 
                             filename = "E:/paper-files/images/go_mf_venn_down_mh_id.tiff",
                             disable.logging = FALSE, 
                             imagetype = "tiff", 
                             main = "Mouse versus Human - Molecular Functions", 
                             main.fontface = "bold", 
                             sub = "Downregulated GO", 
                             col = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                             fill = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                             alpha = 0.3, 
                             fontfamily = "sans", 
                             resolution = 800, 
                             width = 6, 
                             height = 6, 
                             units = "in", 
                             sub.fontface = "bold", 
                             cat.fontfamily = "sans", 
                             cat.fontface = "bold", 
                             cex = 3,
                             cat.cex = 1.5)

#upset plots of the GO down terms
####sets created
#up_bp, down_bp, up_cc, down_cc, up_mf, down_mf
#DOWNREGULATED UPSETPLOTS
CairoPNG("E:/paper-files/images/down_bp_upset.png", width = 6, height = 6, units = "in", res = 800)

upset(
  fromList(down_bp),
  order.by = "freq",  # Order by frequency
  point.size = 3.5,  # Size of points in the matrix
  line.size = 1,  # Size of lines connecting sets
  main.bar.color = "seagreen",  # Color of the main bar
  sets.bar.color = "gray",  # Color of the sets bar
  matrix.color = "royalblue4",  # Color of the matrix points
  text.scale = c(1.8, 1.8, 1.4, 1.4, 1.8, 1.8),  # Scale of text elements
  sets.x.label = "No. of Terms in Set",  # Label for sets bar
  keep.order = TRUE,  # Keep the order of sets
  empty.intersections = "on"  # Show empty intersections
)
# Close the graphics device
dev.off()

CairoPNG("E:/paper-files/images/down_cc_upset.png", width = 6, height = 6, units = "in", res = 800)

upset(
  fromList(down_cc),
  order.by = "freq",  # Order by frequency
  point.size = 3.5,  # Size of points in the matrix
  line.size = 1,  # Size of lines connecting sets
  main.bar.color = "darkorange2",  # Color of the main bar
  sets.bar.color = "gray",  # Color of the sets bar
  matrix.color = "royalblue4",  # Color of the matrix points
  text.scale = c(1.8, 1.8, 1.4, 1.4, 1.8, 1.8),  # Scale of text elements
  sets.x.label = "No. of Terms in Set",  # Label for sets bar
  keep.order = TRUE,  # Keep the order of sets
  empty.intersections = "on"  # Show empty intersections
)
# Close the graphics device
dev.off()

CairoPNG("E:/paper-files/images/down_mf_upset.png", width = 6, height = 6, units = "in", res = 800)

upset(
  fromList(down_mf),
  order.by = "freq",  # Order by frequency
  point.size = 3.5,  # Size of points in the matrix
  line.size = 1,  # Size of lines connecting sets
  main.bar.color = "mediumpurple",  # Color of the main bar
  sets.bar.color = "gray",  # Color of the sets bar
  matrix.color = "royalblue4",  # Color of the matrix points
  text.scale = c(1.8, 1.8, 1.4, 1.4, 1.8, 1.8),  # Scale of text elements
  sets.x.label = "No. of Terms in Set",  # Label for sets bar
  keep.order = TRUE,  # Keep the order of sets
  empty.intersections = "on"  # Show empty intersections
)
# Close the graphics device
dev.off()
##############################TOP TERMS####################################
#filter mRNA by only common GO terms, then filter by most genes and pvalue, qvalue
#need to first make a dataframe of the common elements
# Create a list of vectors
####BIOLOGICAL PROCESSES
vectors_list_up_bp <- list(GO_sxcd_up_venn_bp, GO_sbt_up_venn_bp, GO_adeno_up_venn_bp, GO_hgsoc_up_venn_bp)
elements_up_bp <- calculate.overlap(vectors_list_up_bp)
# Find the maximum length of the vectors
max_length_up_bp <- max(sapply(elements_up_bp, length))
# Pad the vectors with NA to make them of equal length
padded_list_up_bp <- lapply(elements_up_bp, function(vec) {
  length(vec) <- max_length_up_bp
  vec
})
elements_up_df_bp <- data.frame(padded_list_up_bp)
write.csv(elements_up_df_bp, file = "E:/paper-files/venn_up_df_bp.csv", 
          row.names = FALSE)

#repeat for downregulated genes
vectors_list_down_bp <- list(GO_sxcd_down_venn_bp, GO_sbt_down_venn_bp, GO_adeno_down_venn_bp, GO_hgsoc_down_venn_bp)
elements_down_bp <- calculate.overlap(vectors_list_down_bp)
# Find the maximum length of the vectors
max_length_down_bp <- max(sapply(elements_down_bp, length))
# Pad the vectors with NA to make them of equal length
padded_list_down_bp <- lapply(elements_down_bp, function(vec) {
  length(vec) <- max_length_down_bp
  vec
})
elements_down_df_bp <- data.frame(padded_list_down_bp)
write.csv(elements_down_df_bp, file = "E:/paper-files/venn_down_df_bp.csv", 
          row.names = FALSE)

###CELLULAR COMPONENTS
vectors_list_up_cc <- list(GO_sxcd_up_venn_cc, GO_sbt_up_venn_cc, GO_adeno_up_venn_cc, GO_hgsoc_up_venn_cc)
elements_up_cc <- calculate.overlap(vectors_list_up_cc)
# Find the maximum length of the vectors
max_length_up_cc <- max(sapply(elements_up_cc, length))
# Pad the vectors with NA to make them of equal length
padded_list_up_cc <- lapply(elements_up_cc, function(vec) {
  length(vec) <- max_length_up_cc
  vec
})
elements_up_df_cc <- data.frame(padded_list_up_cc)
write.csv(elements_up_df_cc, file = "E:/paper-files/venn_up_df_cc.csv", 
          row.names = FALSE)

#repeat for downregulated genes
vectors_list_down_cc <- list(GO_sxcd_down_venn_cc, GO_sbt_down_venn_cc, GO_adeno_down_venn_cc, GO_hgsoc_down_venn_cc)
elements_down_cc <- calculate.overlap(vectors_list_down_cc)
# Find the maximum length of the vectors
max_length_down_cc <- max(sapply(elements_down_cc, length))
# Pad the vectors with NA to make them of equal length
padded_list_down_cc <- lapply(elements_down_cc, function(vec) {
  length(vec) <- max_length_down_cc
  vec
})
elements_down_df_cc <- data.frame(padded_list_down_cc)
write.csv(elements_down_df_cc, file = "E:/paper-files/venn_down_df_cc.csv", 
          row.names = FALSE)

####-----------------------------------MOLECULAR FUNCTIONS
vectors_list_up_mf <- list(GO_sxcd_up_venn_mf, GO_sbt_up_venn_mf, GO_adeno_up_venn_mf, GO_hgsoc_up_venn_mf)
elements_up_mf <- calculate.overlap(vectors_list_up_mf)
# Find the maximum length of the vectors
max_length_up_mf <- max(sapply(elements_up_mf, length))
# Pad the vectors with NA to make them of equal length
padded_list_up_mf <- lapply(elements_up_mf, function(vec) {
  length(vec) <- max_length_up_mf
  vec
})
elements_up_df_mf <- data.frame(padded_list_up_mf)
write.csv(elements_up_df_mf, file = "E:/paper-files/venn_up_df_mf.csv", 
          row.names = FALSE)

#repeat for downregulated genes
vectors_list_down_mf <- list(GO_sxcd_down_venn_mf, GO_sbt_down_venn_mf, GO_adeno_down_venn_mf, GO_hgsoc_down_venn_mf)
elements_down_mf <- calculate.overlap(vectors_list_down_mf)
# Find the maximum length of the vectors
max_length_down_mf <- max(sapply(elements_down_mf, length))
# Pad the vectors with NA to make them of equal length
padded_list_down_mf <- lapply(elements_down_mf, function(vec) {
  length(vec) <- max_length_down_mf
  vec
})
elements_down_df_mf <- data.frame(padded_list_down_mf)
write.csv(elements_down_df_mf, file = "E:/paper-files/venn_down_df_mf.csv", 
          row.names = FALSE)

####filter GO_mRNA by each common list?
common_go_bp_up  <- na.omit(elements_up_df_bp$a6) #filter for the a6 column which is the intersection
common_go_bp_up <- c(common_go_bp_up) #concatenate into vector
common_go_bp_down  <- na.omit(elements_down_df_bp$a6) #filter for the a6 column which is the intersection
common_go_bp_down <- c(common_go_bp_down)

common_go_cc_up  <- na.omit(elements_up_df_cc$a6) #filter for the a6 column which is the intersection
common_go_cc_up <- c(common_go_cc_up) #concatenate into vector
common_go_cc_down  <- na.omit(elements_down_df_cc$a6) #filter for the a6 column which is the intersection
common_go_cc_down <- c(common_go_cc_down)

common_go_mf_up  <- na.omit(elements_up_df_mf$a6) #filter for the a6 column which is the intersection
common_go_mf_up <- c(common_go_mf_up) #concatenate into vector
common_go_mf_down  <- na.omit(elements_down_df_mf$a6) #filter for the a6 column which is the intersection
common_go_mf_down <- c(common_go_mf_down)

common_go_up <- c(common_go_bp_up, common_go_cc_up, common_go_mf_up)
common_go_down <- c(common_go_bp_down, common_go_cc_down, common_go_mf_down)

common_go_up_df <- GO_mrna %>%
  filter(direction == "Up", ID %in% common_go_up)
common_go_down_df <- GO_mrna %>%
  filter(direction == "Down", ID %in% common_go_down)

common_go <- dplyr::bind_rows(common_go_up_df, common_go_down_df)
write.csv(common_go, file = "E:/paper-files/common_go.csv", row.names = FALSE)

####################Visualise ALL####################
##convert gene ratio into a number i.e. calculate the ratio
#split the character values into two parts
split_result <- strsplit(common_go$GeneRatio, "/")
#Convert the parts to numeric values and perform division
GeneRatioCalc <- sapply(split_result, function(x) as.numeric(x[1]) / as.numeric(x[2]))
#add the numeric result to a new column
common_go$GeneRatioCalc <- GeneRatioCalc

common_go$log10_padj <- -log10(common_go$p.adjust)
write.csv(common_go, "E:/paper-files/common_go.csv", row.names = TRUE)

common_go$group <- factor(common_go$group, levels = c("Sex cords", "Adenoma", "SBT", "HGSOC"))
common_go$direction <- factor(common_go$direction, levels = c("Up", "Down"))
common_go$ONTOLOGY <- factor(common_go$ONTOLOGY, levels = c("BP", "CC", "MF"))

# Arrange the dataframe by the ontology types
common_go <- dplyr::arrange(common_go, ONTOLOGY)

# Set the final factor order for the pathways
unique_descriptions <- unique(common_go$Description)  # Get unique descriptions
common_go$Description <- factor(common_go$Description, levels = rev(unique_descriptions))

gg_common_go <- ggplot(common_go, aes(x = group, y = Description)) + 
  geom_point(aes(fill = log10_padj, size = Count, shape = direction)) + 
  scale_shape_manual(values = c("Up" = 24, "Down" = 25)) + 
  scale_color_gradient(low = "blue", high = "red") + 
  scale_fill_gradient(low = "blue", high = "red") +
  facet_wrap(~ ONTOLOGY) + 
  scale_size_continuous() + 
  labs(x = "", y = "", 
       title = "", 
       subtitle = "",
       color = "'-log10(p.adj)", 
       size = "Gene Count") + 
  theme_grey() +
  theme(
    axis.text.x = element_text(angle = 90, size = 16.0, vjust = 0.5),
    axis.text.y = element_text(size = 16.0, vjust = 0.5),
    axis.title.x = element_text(size = 16.0, vjust = -3.0),
    axis.title.y = element_text(size = 16.0, vjust = 3.0),
    text = element_text(size = 16.0),
    plot.title = element_text(vjust = +3.0, hjust = 0.5),
    plot.margin = margin(1,1,3,1, "cm"),
    legend.position = "bottom",
    legend.margin = margin(t = 1),
    strip.text = element_text(size = 16.0)
  ) +
  guides(
    fill = guide_colorbar(title.position = "top", title.hjust = 0.5),
    size = guide_legend(title.position = "top", title.hjust = 0.5),
    shape = guide_legend(title.position = "top", title.hjust = 0.5)
  )

print(gg_common_go)
ggsave("E:/paper-files/images/common_go.png", plot = gg_common_go, width = 12, height = 14, dpi = 800, 
       bg = "transparent")
ggsave("E:/paper-files/images/common_go_small.png", plot = gg_common_go, width = 12, height = 10, dpi = 300)
#########################################KEGG Analysis##########################
#KEGGenrich requires ENTREZ IDs
#get gene lists again but for ENTREZ ID
sxcd_kegg <- read.csv('E:/paper-files/mlcm_total_sxcd_v_control_deseq_xlfcshrink.csv', sep=',', header = TRUE)
adeno_kegg <- read.csv('E:/paper-files/mlcm_total_adeno_v_control_deseq_xlfcshrink.csv', sep=',', header = TRUE)
sbt_kegg <- read.csv('E:/paper-files/hlcm_total_SBT_v_control_deseq_xlfcshrink.csv', sep=',', header = TRUE)
hgsoc_kegg <- read.csv('E:/paper-files/hlcm_total_HGSOC_v_control_deseq_xlfcshrink.csv', sep=',', header = TRUE)

#filter each dataset for mRNAs
sxcd_kegg_up <- filter(sxcd_kegg, log2FoldChange >= 0, biotype == "protein_coding", entrez != "")
sxcd_kegg_down <- filter(sxcd_kegg, log2FoldChange < 0, biotype == "protein_coding", entrez != "")
adeno_kegg_up <- filter(adeno_kegg, log2FoldChange >= 0, biotype == "protein_coding", entrez != "")
adeno_kegg_down <- filter(adeno_kegg, log2FoldChange < 0, biotype == "protein_coding", entrez != "")
sbt_kegg_up <- filter(sbt_kegg, log2FoldChange >= 0, biotype == "protein_coding", entrez != "")
sbt_kegg_down <- filter(sbt_kegg, log2FoldChange < 0, biotype == "protein_coding", entrez != "")
hgsoc_kegg_up <- filter(hgsoc_kegg, log2FoldChange >= 0, biotype == "protein_coding", entrez != "")
hgsoc_kegg_down <- filter(hgsoc_kegg, log2FoldChange < 0, biotype == "protein_coding", entrez != "")

sxcd_kegg_up <- sxcd_kegg_up$entrez
sxcd_kegg_down <- sxcd_kegg_down$entrez
adeno_kegg_up <- adeno_kegg_up$entrez
adeno_kegg_down <- adeno_kegg_down$entrez
sbt_kegg_up <- sbt_kegg_up$entrez
sbt_kegg_down <- sbt_kegg_down$entrez
hgsoc_kegg_up <- hgsoc_kegg_up$entrez
hgsoc_kegg_down <- hgsoc_kegg_down$entrez
#KEGG ADENO
#make background genes in entrez ID
#make sure the background genes are character vectors
mouse_background <- mouse_background$symbol
human_background <- human_background$symbol

mouse_background_entz <- mapIds(org.Mm.eg.db, keys = mouse_background, 
                                column = "ENTREZID", keytype = "SYMBOL")
mouse_background_entz <- as.character(unlist(mouse_background_entz))
human_background_entz <- mapIds(org.Hs.eg.db, keys = human_background, 
                                column = "ENTREZID", keytype = "SYMBOL")
human_background_entz <- as.character(unlist(human_background_entz))

KEGG_adeno_up <- enrichKEGG(
  gene = adeno_kegg_up, 
  organism = "mmu", 
  keyType = "ncbi-geneid", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  minGSSize = 15, 
  maxGSSize = 500,
  pAdjustMethod = "BH", 
  universe = mouse_background_entz)

KEGG_adeno_up <- as.data.frame(KEGG_adeno_up)
KEGG_adeno_up['direction'] = "Up"

KEGG_adeno_down <- enrichKEGG(
  gene = adeno_kegg_down, 
  organism = "mmu", 
  keyType = "ncbi-geneid", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  minGSSize = 15, 
  maxGSSize = 500,
  pAdjustMethod = "BH", 
  universe = mouse_background_entz)

KEGG_adeno_down <- as.data.frame(KEGG_adeno_down)
KEGG_adeno_down['direction'] = "Down"

#SEXCORD KEGG
KEGG_sxcd_up <- enrichKEGG(
  gene = sxcd_kegg_up, 
  organism = "mmu", 
  keyType = "ncbi-geneid", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  minGSSize = 15, 
  maxGSSize = 500,
  pAdjustMethod = "BH", 
  universe = mouse_background_entz)

KEGG_sxcd_up <- as.data.frame(KEGG_sxcd_up)
KEGG_sxcd_up['direction'] = "Up"

KEGG_sxcd_down <- enrichKEGG(
  gene = sxcd_kegg_down, 
  organism = "mmu", 
  keyType = "ncbi-geneid", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  minGSSize = 15, 
  maxGSSize = 500,
  pAdjustMethod = "BH", 
  universe = mouse_background_entz)

KEGG_sxcd_down <- as.data.frame(KEGG_sxcd_down)
KEGG_sxcd_down['direction'] = "Down"

#SBT KEGG
KEGG_sbt_up <- enrichKEGG(
  gene = sbt_kegg_up, 
  organism = "hsa", 
  keyType = "ncbi-geneid", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  minGSSize = 15, 
  maxGSSize = 500,
  pAdjustMethod = "BH", 
  universe = human_background_entz)

KEGG_sbt_up <- as.data.frame(KEGG_sbt_up)
KEGG_sbt_up['direction'] = "Up"

KEGG_sbt_down <- enrichKEGG(
  gene = sbt_kegg_down, 
  organism = "hsa", 
  keyType = "ncbi-geneid", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  minGSSize = 15, 
  maxGSSize = 500,
  pAdjustMethod = "BH", 
  universe = human_background_entz)

KEGG_sbt_down <- as.data.frame(KEGG_sbt_down)
KEGG_sbt_down['direction'] = "Down"

#HGSOC KEGG
KEGG_hgsoc_up <- enrichKEGG(
  gene = hgsoc_kegg_up, 
  organism = "hsa", 
  keyType = "ncbi-geneid", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  minGSSize = 15, 
  maxGSSize = 500,
  pAdjustMethod = "BH", 
  universe = human_background_entz)

KEGG_hgsoc_up <- as.data.frame(KEGG_hgsoc_up)
KEGG_hgsoc_up['direction'] = "Up"

KEGG_hgsoc_down <- enrichKEGG(
  gene = hgsoc_kegg_down, 
  organism = "hsa", 
  keyType = "ncbi-geneid", 
  pvalueCutoff = 0.05, 
  qvalueCutoff = 0.05, 
  minGSSize = 15, 
  maxGSSize = 500,
  pAdjustMethod = "BH", 
  universe = human_background_entz)

KEGG_hgsoc_down <- as.data.frame(KEGG_hgsoc_down)
KEGG_hgsoc_down['direction'] = "Down"

#combine all GO mouse and human into one data frame and then export
KEGG_adeno_up['group'] = "Adenoma"
KEGG_adeno_up['species'] = "mouse"
KEGG_adeno_down['group'] = "Adenoma"
KEGG_adeno_down['species'] = "mouse"
KEGG_sxcd_up['group'] = "Sex cords"
KEGG_sxcd_up['species'] = "mouse"
KEGG_sxcd_down['group'] = "Sex cords"
KEGG_sxcd_down['species'] = "mouse"
KEGG_sbt_up['group'] = "SBT"
KEGG_sbt_up['species'] = "human"
KEGG_sbt_down['group'] = "SBT"
KEGG_sbt_down['species'] = "human"
KEGG_hgsoc_up['group'] = "HGSOC"
KEGG_hgsoc_up['species'] = "human"
KEGG_hgsoc_down['group'] = "HGSOC"
KEGG_hgsoc_down['species'] = "human"

KEGG_mrna <- dplyr::bind_rows(KEGG_adeno_up, KEGG_adeno_down,
                              KEGG_sxcd_up, KEGG_sxcd_down, 
                              KEGG_sbt_up, KEGG_sbt_down, 
                              KEGG_hgsoc_up, KEGG_hgsoc_down)

write.csv(KEGG_mrna, "E:/paper-files/KEGG_mrna.csv", row.names = TRUE)

#remove mmu and hsa prefix from KEGG IDs so that I can do venn overlaps
KEGG_mrna$ID <- gsub("^(mmu|hsa)", "", KEGG_mrna$ID)
#-----------------------------KEGG VENN--------------
#adenoma KEGG Venn
KEGG_adeno_up_venn <- KEGG_mrna %>%
  dplyr::select(ID, direction, group) %>%
  filter(direction == "Up", group == "Adenoma")
KEGG_adeno_up_venn <- KEGG_adeno_up_venn$ID

KEGG_adeno_down_venn <- KEGG_mrna %>%
  dplyr::select(ID, direction, group) %>%
  filter(direction == "Down", group == "Adenoma")
KEGG_adeno_down_venn <- KEGG_adeno_down_venn$ID

##sxcd KEGG Venn
KEGG_sxcd_up_venn <- KEGG_mrna %>%
  dplyr::select(ID, direction, group) %>%
  filter(direction == "Up", group == "Sex cords")
KEGG_sxcd_up_venn <- KEGG_sxcd_up_venn$ID

KEGG_sxcd_down_venn <- KEGG_mrna %>%
  dplyr::select(ID, direction, group) %>%
  filter(direction == "Down", group == "Sex cords")
KEGG_sxcd_down_venn <- KEGG_sxcd_down_venn$ID

##sbt KEGG Venn
KEGG_sbt_up_venn <- KEGG_mrna %>%
  dplyr::select(ID, direction, group) %>%
  filter(direction == "Up", group == "SBT")
KEGG_sbt_up_venn <- KEGG_sbt_up_venn$ID

KEGG_sbt_down_venn <- KEGG_mrna %>%
  dplyr::select(ID, direction, group) %>%
  filter(direction == "Down", group == "SBT")
KEGG_sbt_down_venn <- KEGG_sbt_down_venn$ID

##hgsoc KEGG Venn
KEGG_hgsoc_up_venn <- KEGG_mrna %>%
  dplyr::select(ID, direction, group) %>%
  filter(direction == "Up", group == "HGSOC")
KEGG_hgsoc_up_venn <- KEGG_hgsoc_up_venn$ID

KEGG_hgsoc_down_venn <- KEGG_mrna %>%
  dplyr::select(ID, direction, group) %>%
  filter(direction == "Down", group == "HGSOC")
KEGG_hgsoc_down_venn <- KEGG_hgsoc_down_venn$ID

#KEGG VENN DIAGRAM VIS
up_kegg <- list(Sexcords = KEGG_sxcd_up_venn, 
              SBT = KEGG_sbt_up_venn, 
              Adenoma = KEGG_adeno_up_venn, 
              HGSOC = KEGG_hgsoc_up_venn)

venn_up_kegg <- venn.diagram(x = up_kegg,
                           category.names = c("Sexcords", "SBT", "Adenoma", "HGSOC"), 
                           filename = "E:/paper-files/images/kegg_venn_up_mh_id.tiff",
                           disable.logging = FALSE, 
                           imagetype = "tiff", 
                           main = "Mouse versus Human - KEGG Pathways", 
                           main.fontface = "bold", 
                           sub = "Upregulated DEGs", 
                           col = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                           fill = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                           alpha = 0.3, 
                           fontfamily = "sans", 
                           resolution = 800, 
                           width = 6, 
                           height = 6, 
                           units = "in", 
                           sub.fontface = "bold", 
                           cat.fontfamily = "sans", 
                           cat.fontface = "bold", 
                           cex = 3,
                           cat.cex = 1.5)

down_kegg <- list(Sexcords = KEGG_sxcd_down_venn, 
                SBT = KEGG_sbt_down_venn, 
                Adenoma = KEGG_adeno_down_venn, 
                HGSOC = KEGG_hgsoc_down_venn)

venn_down_kegg <- venn.diagram(x = down_kegg,
                             category.names = c("Sexcords", "SBT", "Adenoma", "HGSOC"), 
                             filename = "E:/paper-files/images/kegg_venn_down_mh_id.tiff",
                             disable.logging = FALSE, 
                             imagetype = "tiff", 
                             main = "Mouse versus Human - KEGG Pathways", 
                             main.fontface = "bold", 
                             sub = "Downregulated DEGs", 
                             col = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                             fill = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                             alpha = 0.3, 
                             fontfamily = "sans", 
                             resolution = 800, 
                             width = 6, 
                             height = 6, 
                             units = "in", 
                             sub.fontface = "bold", 
                             cat.fontfamily = "sans", 
                             cat.fontface = "bold", 
                             cex = 3,
                             cat.cex = 1.5)

CairoPNG("E:/paper-files/images/up_kegg_upset.png", width = 6, height = 6, units = "in", res = 800)

upset(
  fromList(up_kegg),
  order.by = "freq",  # Order by frequency
  point.size = 3.5,  # Size of points in the matrix
  line.size = 1,  # Size of lines connecting sets
  main.bar.color = "turquoise4",  # Color of the main bar
  sets.bar.color = "gray",  # Color of the sets bar
  matrix.color = "red3",  # Color of the matrix points
  text.scale = c(1.8, 1.8, 1.4, 1.4, 1.8, 1.8),  # Scale of text elements
  sets.x.label = "No. of Pathways in Set",  # Label for sets bar
  keep.order = TRUE,  # Keep the order of sets
  empty.intersections = "on"  # Show empty intersections
)
# Close the graphics device
dev.off()

CairoPNG("E:/paper-files/images/down_kegg_upset.png", width = 6, height = 6, units = "in", res = 800)

upset(
  fromList(down_kegg),
  order.by = "freq",  # Order by frequency
  point.size = 3.5,  # Size of points in the matrix
  line.size = 1,  # Size of lines connecting sets
  main.bar.color = "turquoise4",  # Color of the main bar
  sets.bar.color = "gray",  # Color of the sets bar
  matrix.color = "royalblue4",  # Color of the matrix points
  text.scale = c(1.8, 1.8, 1.4, 1.4, 1.8, 1.8),  # Scale of text elements
  sets.x.label = "No. of Pathways in Set",  # Label for sets bar
  keep.order = TRUE,  # Keep the order of sets
  empty.intersections = "on"  # Show empty intersections
)
# Close the graphics device
dev.off()

################################KEGG UP VENN ELEMENTS###########################
vectors_list_up_kegg <- list(KEGG_sxcd_up_venn, KEGG_adeno_up_venn, KEGG_sbt_up_venn, KEGG_hgsoc_up_venn)
elements_up_kegg <- calculate.overlap(vectors_list_up_kegg)
# Find the maximum length of the vectors
max_length_up_kegg <- max(sapply(elements_up_kegg, length))
# Pad the vectors with NA to make them of equal length
padded_list_up_kegg <- lapply(elements_up_kegg, function(vec) {
  length(vec) <- max_length_up_kegg
  vec
})
elements_up_df_kegg <- data.frame(padded_list_up_kegg)
write.csv(elements_up_df_kegg, file = "E:/paper-files/venn_up_df_kegg.csv", 
          row.names = FALSE)

#repeat for downregulated genes
vectors_list_down_kegg <- list(KEGG_sxcd_down_venn, KEGG_adeno_down_venn, KEGG_sbt_down_venn, KEGG_hgsoc_down_venn)
elements_down_kegg <- calculate.overlap(vectors_list_down_kegg)
# Find the maximum length of the vectors
max_length_down_kegg <- max(sapply(elements_down_kegg, length))
# Pad the vectors with NA to make them of equal length
padded_list_down_kegg <- lapply(elements_down_kegg, function(vec) {
  length(vec) <- max_length_down_kegg
  vec
})
elements_down_df_kegg <- data.frame(padded_list_down_kegg)
write.csv(elements_down_df_kegg, file = "E:/paper-files/venn_down_df_kegg.csv", 
          row.names = FALSE)

#-------------------Plot KEGG
#there are no common and not that many terms soooo plot all??
KEGG_mrna$Description <- gsub(" - Mus musculus \\(house mouse\\)", "", KEGG_mrna$Description)


#######################PLOT KEGG - TOP########################################
top_10_kegg <- KEGG_mrna %>%
  arrange(group, direction, desc(Count)) %>%
  group_by(group, direction) %>%
  slice_head(n = 10)

write.csv(top_10_kegg, file = "E:/paper-files/top_10_kegg.csv", 
          row.names = FALSE)

top_5_kegg <- KEGG_mrna %>%
  arrange(group, direction, desc(Count)) %>%
  group_by(group, direction) %>%
  slice_head(n = 5)

top_10_kegg <- top_10_kegg %>%
  arrange(subcategory, Description) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

write.csv(top_5_kegg, file = "E:/paper-files/top_5_kegg.csv", 
          row.names = FALSE)

##convert gene ratio into a number i.e. calculate the ratio
#split the character values into two parts
split_result <- strsplit(top_10_kegg$GeneRatio, "/")
#Convert the parts to numeric values and perform division
GeneRatioCalc <- sapply(split_result, function(x) as.numeric(x[1]) / as.numeric(x[2]))
#add the numeric result to a new column
top_10_kegg$GeneRatioCalc <- GeneRatioCalc

top_10_kegg$log10_padj <- -log10(top_10_kegg$p.adjust)
write.csv(top_10_kegg, "E:/paper-files/top_10_kegg.csv", row.names = TRUE)

top_10_kegg$group <- factor(top_10_kegg$group, levels = c("sex-cords", "adenoma", "SBT", "HGSOC"))
top_10_kegg$direction <- factor(top_10_kegg$direction, levels = c("Up", "Down"))

# Arrange the dataframe by the ontology types
#common_go <- dplyr::arrange(common_go, ONTOLOGY)

# Set the final factor order for the pathways

#unique_descriptions <- unique(common_go$Description)  # Get unique descriptions
#common_go$Description <- factor(common_go$Description, levels = rev(unique_descriptions))

# Wrap long facet labels
wrapped_labeller <- labeller(subcategory = label_wrap_gen(width = 25))

gg_top_10_kegg <- ggplot(top_10_kegg, aes(x = group, y = Description)) + 
  geom_point(aes(size = Count, shape = direction, fill = log10_padj)) + 
  facet_grid(subcategory ~ ., scales = "free_y", space = "free_y", labeller = wrapped_labeller) + 
  scale_size_continuous() + 
  scale_color_gradient(low = "blue", high = "red") + 
  scale_fill_gradient(low = "blue", high = "red") +
  scale_shape_manual(values = c("Up" = 24, "Down" = 25)) + 
  labs(x = "Group", y = "KEGG Pathway", 
       title = "KEGG Pathway Analysis", 
       subtitle = "for Significant DEGs (Group/Control)",
       fill = "-log10(p.adj)", 
       size = "Gene Count") + 
  theme_grey() +
  theme(
    axis.text.x = element_text(angle = 90, size = 18.0, vjust = 0.5),
    axis.text.y = element_text(size = 18.0, vjust = 0.5),
    axis.title.x = element_text(size = 18.0, vjust = -3.0),
    axis.title.y = element_text(size = 18.0, vjust = 3.0),
    text = element_text(size = 18.0),
    plot.title = element_text(vjust = +3.0, hjust = 0.5),
    plot.margin = margin(1,1,1,1, "cm"), 
    strip.text.y = element_text(size = 17, angle = 0, hjust = 0.5)
  )

print(gg_top_10_kegg)
ggsave("E:/paper-files/images/top_10_kegg_small.png", plot = gg_top_10_kegg, width = 16, height = 19, dpi = 800)

#---------------------------KEGG overlap irrespective of Up or Down direction
####-----------------------------------KEGG UP VENN ELEMENTS
duplicates <- KEGG_mrna %>%
  group_by(ID) %>%
  filter(n() > 1) %>%
  distinct(ID, .keep_all = TRUE)

duplicates <- duplicates$ID

KEGG_mrna_dups <- KEGG_mrna %>%
  filter(ID %in% duplicates)

write.csv(KEGG_mrna_dups, file = "E:/paper-files/KEGG_mrna_dups.csv", 
          row.names = FALSE)

##convert gene ratio into a number i.e. calculate the ratio
#split the character values into two parts
split_result <- strsplit(KEGG_mrna_dups$GeneRatio, "/")
#Convert the parts to numeric values and perform division
GeneRatioCalc <- sapply(split_result, function(x) as.numeric(x[1]) / as.numeric(x[2]))
#add the numeric result to a new column
KEGG_mrna_dups$GeneRatioCalc <- GeneRatioCalc

KEGG_mrna_dups$log10_padj <- -log10(KEGG_mrna_dups$p.adjust)
write.csv(KEGG_mrna_dups, "E:/paper-files/KEGG_mrna_dups.csv", row.names = TRUE)

KEGG_mrna_dups$group <- factor(KEGG_mrna_dups$group, levels = c("Sex cords", "Adenoma", "SBT", "HGSOC"))
KEGG_mrna_dups$direction <- factor(KEGG_mrna_dups$direction, levels = c("Up", "Down"))

# Wrap long facet labels
wrapped_labeller <- labeller(subcategory = label_wrap_gen(width = 25))

kegg_dups <- ggplot(KEGG_mrna_dups, aes(x = group, y = Description)) + 
  geom_point(aes(size = Count, shape = direction, fill = log10_padj)) + 
  facet_grid(subcategory ~ ., scales = "free_y", space = "free_y",  labeller = label_wrap_gen(width = 35)) + 
  scale_size_continuous() + 
  scale_color_gradient(low = "blue", high = "red") + 
  scale_fill_gradient(low = "blue", high = "red") +
  scale_shape_manual(values = c("Up" = 24, "Down" = 25)) + 
  labs(x = "", y = "", 
       title = "KEGG Pathway Analysis", 
       subtitle = "for Significant DEGs (Group/Control)",
       fill = "-log10(p.adj)", 
       size = "Gene Count") + 
  theme_grey() +
  theme(
    axis.text.x = element_text(angle = 90, size = 20.0, vjust = 0.5),
    axis.text.y = element_text(size = 18.0, vjust = 0.5),
    axis.title.x = element_text(size = 20.0, vjust = -3.0),
    axis.title.y = element_text(size = 20.0, vjust = 3.0),
    text = element_text(size = 16.0),
    plot.title = element_text(vjust = +3.0, hjust = 0.5),
    plot.margin = margin(1,1,1,1, "cm"), 
    strip.text.y = element_text(size = 17, angle = 0, hjust = 0.5), 
    panel.spacing.y = unit(0.2, "lines"),  
   legend.position = "right",  legend.margin = margin(t = 5),
  )

print(kegg_dups)
ggsave("E:/paper-files/images/KEGG_mrna_dups_test.png", plot = kegg_dups, width = 16, height = 20, 
       dpi = 800, bg = "transparent")

#----------------------------------GO ENRICH COMMON DEGS------------------------
#refresh dataset for correct symbol types
sxcd <- read.csv('E:/paper-files/mlcm_total_sxcd_v_control_deseq_xlfcshrink.csv', sep=',', header = TRUE)
sbt <- read.csv('E:/paper-files/hlcm_total_SBT_v_control_deseq_xlfcshrink.csv', sep=',', header = TRUE)

#get common gene lists
common_degs_up  <- na.omit(elements_up_df$a6) #filter for the a6 column which is the intersection
common_degs_up <- c(common_degs_up) #concatenate into vector
common_degs_down  <- na.omit(elements_down_df$a6) #filter for the a6 column which is the intersection
common_degs_down <- c(common_degs_down) #concatenate into vector

#get correct symbol notation 'SYMBOL' by filtering DEG dataframe
#end up with lists of Gene case for  mouse and GENE case for human
mouse_common_degs_up <- sxcd %>%
  filter(casematch %in% common_degs_up)
mouse_common_degs_up <- mouse_common_degs_up$symbol
mouse_common_degs_up <- c(mouse_common_degs_up)

mouse_common_degs_down <- sxcd %>%
  filter(casematch %in% common_degs_down)
mouse_common_degs_down <- mouse_common_degs_down$symbol
mouse_common_degs_down <- c(mouse_common_degs_down)

human_common_degs_up <- sbt %>%
  filter(casematch %in% common_degs_up)
human_common_degs_up <- human_common_degs_up$symbol
human_common_degs_up <- c(human_common_degs_up)

human_common_degs_down <- sbt %>%
  filter(casematch %in% common_degs_down)
human_common_degs_down <- human_common_degs_down$symbol
human_common_degs_down <- c(human_common_degs_down)

#DO GOenrich for mouse common genes, with the same background lists
GO_mouse_up <- enrichGO(gene = mouse_common_degs_up, OrgDb = "org.Mm.eg.db", 
                        keyType = "SYMBOL", ont = "ALL", 
                        pAdjustMethod = "fdr", 
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.05, 
                        minGSSize = 15, 
                        maxGSSize = 500, 
                        universe = mouse_background)

GO_mouse_down <- enrichGO(gene = mouse_common_degs_down, OrgDb = "org.Mm.eg.db", 
                        keyType = "SYMBOL", ont = "ALL", 
                        pAdjustMethod = "fdr", 
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.05, 
                        minGSSize = 15, 
                        maxGSSize = 500, 
                        universe = mouse_background)

GO_human_up <- enrichGO(gene = human_common_degs_up, OrgDb = "org.Hs.eg.db", 
                        keyType = "SYMBOL", ont = "ALL", 
                        pAdjustMethod = "fdr", 
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.05, 
                        minGSSize = 15, 
                        maxGSSize = 500, 
                        universe = human_background)

GO_human_down <- enrichGO(gene = human_common_degs_down, OrgDb = "org.Hs.eg.db", 
                          keyType = "SYMBOL", ont = "ALL", 
                          pAdjustMethod = "fdr", 
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.05, 
                          minGSSize = 15, 
                          maxGSSize = 500, 
                          universe = human_background)

GO_mouse_up <- as.data.frame(GO_mouse_up)
GO_mouse_up['direction'] = "Up"
GO_mouse_up['species'] = "mouse"

GO_mouse_down <- as.data.frame(GO_mouse_down)
GO_mouse_down['direction'] = "Down"
GO_mouse_down['species'] = "mouse"

GO_human_up <- as.data.frame(GO_human_up)
GO_human_up['direction'] = "Up"
GO_human_up['species'] = "human"

GO_human_down <- as.data.frame(GO_human_down)
GO_human_down['direction'] = "Down"
GO_human_down['species'] = "human"

GO_common_degs <- dplyr::bind_rows(GO_mouse_up, GO_mouse_down, GO_human_up, GO_human_down)
write.csv(GO_common_degs, "E:/paper-files/GO_common_genes.csv", row.names = TRUE)

##convert gene ratio into a number i.e. calculate the ratio
#split the character values into two parts
split_result <- strsplit(GO_common_degs$GeneRatio, "/")
#Convert the parts to numeric values and perform division
GeneRatioCalc <- sapply(split_result, function(x) as.numeric(x[1]) / as.numeric(x[2]))
#add the numeric result to a new column
GO_common_degs$GeneRatioCalc <- GeneRatioCalc

GO_common_degs$log10_padj <- -log10(GO_common_degs$p.adjust)
write.csv(GO_common_degs, "E:/paper-files/GO_common_degs.csv", row.names = TRUE)

GO_common_degs$direction <- factor(GO_common_degs$direction, levels = c("Up", "Down"))
GO_common_degs$ONTOLOGY <- factor(GO_common_degs$ONTOLOGY, levels = c("BP", "CC", "MF"))
GO_common_degs$species <- factor(GO_common_degs$species, levels = c("mouse", "human"))

# Arrange the dataframe by the ontology types
GO_common_degs <- dplyr::arrange(GO_common_degs, ONTOLOGY)

# Set the final factor order for the pathways
unique_descriptions <- unique(GO_common_degs$Description)  # Get unique descriptions
GO_common_degs$Description <- factor(GO_common_degs$Description, levels = rev(unique_descriptions))

#filter out terms that have 10 or more genes associated:
GO_common_degs_10 <- GO_common_degs %>%
  filter(Count >= 10)

GO_common_degs_5 <- GO_common_degs %>%
  filter(Count >= 5)

GO_common_degs_down <- GO_common_degs %>%
  filter(direction == "Down")


gg_GO_common_degs <- ggplot(GO_common_degs, aes(x = species, y = Description)) + 
  geom_point(aes(color = log10_padj, size = Count, shape = direction)) + 
  facet_wrap(~ ONTOLOGY) + 
  scale_size_continuous() + 
  scale_color_gradient(low = "blue", high = "red") + 
  labs(x = "Group", y = "Gene Ontology Term", 
       title = "Gene Ontology Analysis", 
       subtitle = "for Significant DEGs (Group/Control)",
       color = "-log10(p.adj)", 
       size = "Gene Count") + 
  theme_gray() +
  theme(
    axis.text.x = element_text(angle = 90, size = 12.0, vjust = 0.5),
    axis.text.y = element_text(size = 12.0, vjust = 0.5),
    axis.title.x = element_text(size = 13.0, vjust = -3.0),
    axis.title.y = element_text(size = 13.0, vjust = 3.0),
    text = element_text(size = 14.0),
    plot.title = element_text(vjust = +3.0, hjust = 0.5),
    plot.margin = margin(1,1,1,1, "cm")
  )

print(gg_GO_common_degs)

ggsave("E:/paper-files/images/GO_common_degs.png", plot = gg_GO_common_degs, 
       width = 12, height = 24, dpi = 800)
ggsave("E:/paper-files/images/GO_common_degs_small.png", plot = gg_GO_common_degs, 
       width = 12, height = 10, dpi = 300)
#---------------------repeat for only top 10------------------------------------
gg_GO_common_degs_5 <- ggplot(GO_common_degs_5, aes(x = species, y = Description)) + 
  geom_point(aes(color = log10_padj, size = Count, shape = direction)) + 
  facet_wrap(~ ONTOLOGY) + 
  scale_size_continuous() + 
  scale_color_gradient(low = "blue", high = "red") + 
  labs(x = "Group", y = "Gene Ontology Term", 
       title = "Gene Ontology Analysis", 
       subtitle = "for Significant DEGs (Group/Control)",
       color = "-log10(p.adj)", 
       size = "Gene Count") + 
  theme_gray() +
  theme(
    axis.text.x = element_text(angle = 90, size = 12.0, vjust = 0.5),
    axis.text.y = element_text(size = 12.0, vjust = 0.5),
    axis.title.x = element_text(size = 13.0, vjust = -3.0),
    axis.title.y = element_text(size = 13.0, vjust = 3.0),
    text = element_text(size = 14.0),
    plot.title = element_text(vjust = +3.0, hjust = 0.5),
    plot.margin = margin(1,1,1,1, "cm")
  )

print(gg_GO_common_degs_5)

ggsave("E:/paper-files/images/GO_common_degs_5.png", plot = gg_GO_common_degs_5, 
       width = 12, height = 20, dpi = 800)
ggsave("E:/paper-files/images/GO_common_degs_5_small.png", plot = gg_GO_common_degs_5, 
       width = 10, height = 12, dpi = 300)

gg_GO_common_degs_down <- ggplot(GO_common_degs_down, aes(x = species, y = Description)) + 
  geom_point(aes(color = log10_padj, size = Count, shape = direction)) + 
  facet_wrap(~ ONTOLOGY) + 
  scale_size_continuous(breaks = c(3, 4)) + 
  scale_color_gradient(low = "blue", high = "red") + 
  labs(x = "Group", y = "Gene Ontology Term", 
       title = "Gene Ontology Analysis", 
       subtitle = "for Common Significant DEGs",
       color = "-log10(p.adj)", 
       size = "Gene Count") + 
  theme_gray() +
  theme(
    axis.text.x = element_text(angle = 90, size = 12.0, vjust = 0.5),
    axis.text.y = element_text(size = 12.0, vjust = 0.5),
    axis.title.x = element_text(size = 13.0, vjust = -3.0),
    axis.title.y = element_text(size = 13.0, vjust = 3.0),
    text = element_text(size = 14.0),
    plot.title = element_text(vjust = +3.0, hjust = 0.5),
    plot.margin = margin(1,1,1,1, "cm")
  )

print(gg_GO_common_degs_down)

ggsave("E:/paper-files/images/GO_common_degs_down.png", plot = gg_GO_common_degs_down, 
       width = 6, height = 5, dpi = 800)
ggsave("E:/paper-files/images/GO_common_degs_down_small.png", plot = gg_GO_common_degs_down, 
       width = 10, height = 12, dpi = 300)
