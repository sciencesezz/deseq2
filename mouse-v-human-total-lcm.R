
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


#############################Upload files as variables##########################
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

sxcd_up <- arrange(sxcd_up, -log2FoldChange) %>%
  select(symbol)
sxcd_down <- arrange(sxcd_down, log2FoldChange) %>%
  select(symbol)
adeno_up <- arrange(adeno_up, -log2FoldChange) %>%
  select(symbol)
adeno_down <- arrange(adeno_down, log2FoldChange) %>%
  select(symbol)
sbt_up <- arrange(sbt_up, -log2FoldChange) %>%
  select(symbol)
sbt_down <- arrange(sbt_down, log2FoldChange) %>%
  select(symbol)
hgsoc_up <- arrange(hgsoc_up, -log2FoldChange) %>%
  select(symbol)
hgsoc_down <- arrange(hgsoc_down, log2FoldChange) %>%
  select(symbol)

#coerce to vector
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
                             filename = "E:/paper-files/images/venn_up_mh.tiff",
                             disable.logging = FALSE, 
                             imagetype = "tiff", 
                             main = "Mouse versus Human", 
                             main.fontface = "bold", 
                             sub = "Upregulated DEGs", 
                             col = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                             fill = c("palegreen1", "mediumpurple1","hotpink", "grey"), 
                             alpha = 0.3, 
                          fontfamily = "sans", 
                          resolution = 600, 
                          width = 6, 
                          height = 6, 
                          units = "in", 
                          sub.fontface = "bold", 
                          cat.fontfamily = "sans", 
                          cat.fontface = "bold", 
                          cex = 2,
                          cat.cex = 1.5)

venn_down_mh <- venn.diagram(x = venn_data_ordered_down,
                           category.names = c("Sexcords", "SBT", "Adenoma", "HGSOC"), 
                           filename = "E:/paper-files/images/venn_down_mh.tiff",
                           disable.logging = FALSE, 
                           imagetype = "tiff", 
                           main = "Mouse versus Human", 
                           main.fontface = "bold", 
                           sub = "Downregulated DEGs", 
                           col = c("palegreen1", "mediumpurple1", "hotpink", "grey"), 
                           fill = c("palegreen1", "mediumpurple1","hotpink", "grey"), 
                           alpha = 0.3, 
                           fontfamily = "sans", 
                           resolution = 600, 
                           width = 6, 
                           height = 6, 
                           units = "in", 
                           sub.fontface = "bold", 
                           cat.fontfamily = "sans", 
                           cat.fontface = "bold", 
                           cex = 2, 
                           cat.cex = 1.5)

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
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE) {
  gmt <- gmtPathways(gmt_file)
  hidden <- unique(unlist(gmt))
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)), 
                nrow = length(hidden), ncol = length(gmt))
  for (i in 1:dim(mat)[2]){
    mat[,i] <- as.numeric(hidden %in% gmt[[i]])
  }
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1,])>5)]]
  final_list <- matrix_to_list(mat)
  if(savefile){
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('Wohoo! .gmt conversion successfull!:)')
  return(final_list)
}


sxcd_go <- as.data.frame(sxcd_lfc)
sxcd_go <- sxcd_go %>%
  mutate(symbol = str_to_title(symbol))
sxcd_go$log2FoldChange <- as.numeric(sxcd_go$log2FoldChange)
sxcd_go$padj <- as.numeric(sxcd_go$padj)
sxcd_go <- sxcd_go%>%
  arrange(dplyr::desc(log2FoldChange))
rankings <- sign(sxcd_go$log2FoldChange)*(-log10(sxcd_go$padj))
names(rankings) <- sxcd_go$symbol


go_pathway <- gmtPathways("E:/paper-files/m5.go.v2023.2.Mm.symbols.gmt")
tumour_pathway <- gmtPathways("E:/paper-files/m5.mpt.v2023.2.Mm.symbols.gmt")
cc_pathway <- gmtPathways("E:/paper-files/m5.go.cc.v2023.2.Mm.symbols.gmt")


mygenes_mouse <- read.csv('E:/paper-files/mlcm_total_logcpmc.csv', sep=',', header = TRUE)
colnames(mygenes_mouse)[1] <- "ensgene"
mygenes_mouse <- left_join(x = mygenes_mouse,
                          y = grcm38 [, c("ensgene", "symbol", "entrez", "biotype", "description")], 
                          by = "ensgene")
mygenes_mouse <- subset(mygenes_mouse, biotype == "protein_coding")
mygenes_mouse <- mygenes_mouse$symbol


# Run fgsea
fgsea_sxcd_go <- fgsea(pathways = go_pathway, 
                       stats = rankings ,
                       minSize = min_size,
                       maxSize = max_size)

fgsea_sxcd_cc <- fgseaMultilevel(pathways = cc_pathway, 
                  stats = rankings ,
                  minSize = min_size,
                  maxSize = max_size)

fgsea_sxcd_tumour <- fgsea(pathways = tumour_pathway, 
                    stats = sxcd_go ,
                    minSize = min_size,
                    maxSize = max_size)


GSEA_sxcd <- gseGO(geneList = geneList, 
                   ont = "ALL", 
                   OrgDb = "org.Mm.eg.db", 
                   keyType = "SYMBOL", 
                   minGSSize = min_size, 
                   maxGSSize = max_size, 
                   pvalueCutoff = padj_cutoff, 
                   verbose = TRUE) 

GSEA_sxcd <- as.data.frame(GSEA_sxcd)
GSEA_adeno <- as.data.frame(GSEA_adeno)
GSEA_sbt <- as.data.frame(GSEA_sbt)
GSEA_hgsoc <- as.data.frame(GSEA_hgsoc)

GSEA_sxcd['sample'] <- "sex-cords"
GSEA_adeno['sample'] <- "adenoma"
GSEA_sbt['sample'] <- "SBT"
GSEA_hgsoc['sample'] <- "HGSOC"

GSEA_all<- dplyr::bind_rows(GSEA_sxcd, GSEA_adeno, GSEA_sbt, GSEA_hgsoc)

write.csv(GSEA_all, "E:/GSEA_all.csv", row.names = TRUE)

#filter dataframe to have descriptions that are present in all groups
GSEA_all_filtered <- GSEA_all %>%
  group_by(Description) %>%
  filter(n() == 4) %>%
  ungroup()

GSEA_all_filtered$log10_padj <- -log10(GSEA_all_filtered$p.adjust)

write.csv(GSEA_all_filtered, "E:/GSEA_all_filtered.csv", row.names = TRUE)

GSEA_all_filtered$sample <- factor(GSEA_all_filtered$sample, levels = c("sex-cords", "adenoma", "SBT", "HGSOC"))

GSEA_all_filtered <- GSEA_all_filtered %>% 
    arrange(ONTOLOGY)

GSEA_all_filtered_graph <- ggplot(GSEA_all_filtered, aes(x = sample, y = Description)) + 
  geom_point(aes(color = NES, size = log10_padj, shape = ONTOLOGY)) + 
  #facet_wrap(~ ONTOLOGY, scales = "free_y", ncol = 1, space = "free_y") + 
    scale_size_continuous() + 
  scale_color_gradient(low = "blue", high = "red") + 
  labs(x = "Group", y = "GO Term", 
       title = "Gene Set Enrichment Analysis", 
       subtitle = "Based on Significantly Differentially Genes (Group vs. Control)",
       color = "NES", 
       size = "-log10(p.adj)", 
       shape = "Ontology") + 
  theme_gray() +
  theme(
    axis.text.x = element_text(angle = 0, size = 12.0, vjust = 0.5),
    axis.text.y = element_text(size = 12.0, vjust = 0.5),
    axis.title.x = element_text(size = 13.0, vjust = -3.0),
    axis.title.y = element_text(size = 13.0, vjust = 3.0),
    text = element_text(size = 12.0),
    plot.title = element_text(vjust = +3.0, hjust = 0.5),
    plot.subtitle = element_text(vjust = +3.0, hjust = 0.5), 
    plot.margin = margin(1,1,1,1, "cm")
  )


print(GSEA_all_filtered_graph)


ggsave("E:/GSEA_all_filtered_graph_3.png", plot = GSEA_all_filtered_graph, width = 12, height = 12, dpi = 600)

top_GSEA_sxcd <- GSEA_sxcd %>%
  arrange(ONTOLOGY, p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 5)

top_GSEA_adeno <- GSEA_adeno %>%
  arrange(ONTOLOGY, p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 5)


top_GSEA_sbt <- GSEA_sbt %>%
  arrange(ONTOLOGY, p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 5)


top_GSEA_hgsoc <- GSEA_hgsoc %>%
  arrange(ONTOLOGY, p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 5)


GSEA_all_top <- dplyr::bind_rows(top_GSEA_sxcd, top_GSEA_adeno, top_GSEA_sbt, top_GSEA_hgsoc)

GSEA_all$log10_padj <- -log10(GSEA_all$p.adjust)
GSEA_all$log10_padj <- -log10(GSEA_all$p.adjust)

write.csv(GSEA_all, "E:/GSEA_all_top_5.csv", row.names = TRUE)

GSEA_all <- GSEA_all %>% 
  arrange(NES)
GSEA_all <- subset(GSEA_all, ONTOLOGY != "NA")

GSEA_all$sample <- factor(GSEA_all$sample, levels = c("sex-cords", "adenoma", "SBT", "HGSOC"))


print(GSEA_all)



GSEA_all_graph <- ggplot(GSEA_all, aes(x = NES, y = Description)) + 
  geom_point(aes(color = log10_padj, size = setSize, shape = sample)) + 
  facet_wrap(~ONTOLOGY) + 
  scale_size_continuous() + 
  scale_color_gradient(low = "green", high = "purple") + 
  labs(x = "Group", y = "GO Term", 
       title = "Gene Ontology Enrichment Analysis", 
       subtitle = "Significantly Differentially Upregulated Genes (Group/Control)",
       color = "-log10(p.adj)", 
       size = "Gene Count", 
       shape = "Sample Group") + 
  theme_gray() +
  theme(
    axis.text.x = element_text(angle = 45, size = 10.0, vjust = 0.5),
    axis.text.y = element_text(size = 10.0, vjust = 0.5),
    axis.title.x = element_text(size = 13.0, vjust = -3.0),
    axis.title.y = element_text(size = 13.0, vjust = 3.0),
    text = element_text(size = 10.0),
    plot.title = element_text(vjust = +3.0, hjust = 0.5),
    plot.margin = margin(1,1,1,1, "cm")
  )

print(GSEA_all_graph)

ggsave("E:/GSEA_all_graph.png", plot = GSEA_all_graph, width = 12, height = 12, dpi = 600)

ridgeplot(GSEA_all)

###############################ADENO###########################################


GO_adeno_up <- enrichGO(gene = adeno_up, OrgDb = "org.Mm.eg.db", 
                       keyType = "SYMBOL", ont = "ALL", 
                       pAdjustMethod = "fdr", 
                       pvalueCutoff = padj_cutoff, 
                       qvalueCutoff = qval_cutoff, 
                       minGSSize = min_size, 
                       maxGSSize = max_size)
GO_adeno_up <- as.data.frame(GO_adeno_up)
GO_adeno_up['direction'] = "Up"

GO_adeno_down <- enrichGO(gene = adeno_down, OrgDb = "org.Mm.eg.db", 
                         keyType = "SYMBOL", ont = "ALL", 
                         pAdjustMethod = "fdr", 
                         pvalueCutoff = padj_cutoff, 
                         qvalueCutoff = qval_cutoff, 
                         minGSSize = min_size, 
                         maxGSSize = max_size)
GO_adeno_down <- as.data.frame(GO_adeno_down)
GO_adeno_down['direction'] = "Down"

GO_adeno <- dplyr::bind_rows(GO_adeno_up, GO_adeno_down)
write.csv(GO_adeno, "E:/GO_adeno.csv", row.names = TRUE)
########################################HUMAN GO################################
GO_sbt_up <- enrichGO(gene = sbt_up, OrgDb = "org.Hs.eg.db", 
                            keyType = "SYMBOL", ont = "ALL", 
                            pAdjustMethod = "fdr", 
                            pvalueCutoff = padj_cutoff, 
                            qvalueCutoff = qval_cutoff, 
                            minGSSize = min_size, 
                            maxGSSize = max_size)
GO_sbt_up <- as.data.frame(GO_sbt_up)
GO_sbt_up['direction'] = "Up"

GO_sbt_down <- enrichGO(gene = sbt_down, OrgDb = "org.Hs.eg.db", 
                       keyType = "SYMBOL", ont = "ALL", 
                       pAdjustMethod = "fdr", 
                       pvalueCutoff = padj_cutoff, 
                       qvalueCutoff = qval_cutoff, 
                       minGSSize = min_size, 
                       maxGSSize = max_size)
GO_sbt_down <- as.data.frame(GO_sbt_down)
GO_sbt_down['direction'] = "Down"

GO_sbt <- dplyr::bind_rows(GO_sbt_up, GO_sbt_down)
write.csv(GO_sxcd, "E:/GO_sbt.csv", row.names = TRUE)

GO_hgsoc_up <- enrichGO(gene = hgsoc_up, OrgDb = "org.Hs.eg.db", 
                       keyType = "SYMBOL", ont = "ALL", 
                       pAdjustMethod = "fdr", 
                       pvalueCutoff = padj_cutoff, 
                       qvalueCutoff = qval_cutoff, 
                       minGSSize = min_size, 
                       maxGSSize = max_size)
GO_hgsoc_up <- as.data.frame(GO_hgsoc_up)
GO_hgsoc_up['direction'] = "Up"

GO_hgsoc_down <- enrichGO(gene = hgsoc_down, OrgDb = "org.Hs.eg.db", 
                         keyType = "SYMBOL", ont = "ALL", 
                         pAdjustMethod = "fdr", 
                         pvalueCutoff = padj_cutoff, 
                         qvalueCutoff = qval_cutoff, 
                         minGSSize = min_size, 
                         maxGSSize = max_size)
GO_hgsoc_down <- as.data.frame(GO_hgsoc_down)
GO_hgsoc_down['direction'] = "Down"

GO_hgsoc <- dplyr::bind_rows(GO_hgsoc_up, GO_hgsoc_down)
write.csv(GO_hgsoc, "E:/GO_hgsoc.csv", row.names = TRUE)
#combine all GO mouse and human into one data frame and then export

GO_sxcd['group'] = "sex-cords"
GO_sxcd['species'] = "mouse"
GO_adeno['group'] = "adenoma"
GO_adeno['species'] = "mouse"
GO_sbt['group'] = "SBT"
GO_sbt['species'] = "human"
GO_hgsoc['group'] = "HGSOC"
GO_hgsoc['species'] = "human"

GO_mrna <- dplyr::bind_rows(GO_sxcd, GO_adeno, GO_sbt, GO_hgsoc)
write.csv(GO_mrna, "E:/GO_mrna.csv", row.names = TRUE)
#plotting top 5 terms of each
top_terms_sxcd <- GO_sxcd %>%
  arrange(direction, ONTOLOGY, p.adjust) %>%
  group_by(ONTOLOGY, direction) %>%
  slice_head(n = 5)

top_terms_adeno <- GO_adeno %>%
  arrange(direction, ONTOLOGY, p.adjust) %>%
  group_by(ONTOLOGY, direction) %>%
  slice_head(n = 5)

top_terms_sbt <- GO_sbt %>%
  arrange(direction, ONTOLOGY, p.adjust) %>%
  group_by(ONTOLOGY, direction) %>%
  slice_head(n = 5)

top_terms_hgsoc <- GO_hgsoc %>%
  arrange(direction, ONTOLOGY, p.adjust) %>%
  group_by(ONTOLOGY, direction) %>%
  slice_head(n = 5)

GO_mrna_top <- dplyr::bind_rows(top_terms_sxcd, top_terms_adeno, top_terms_sbt, top_terms_hgsoc)

##convert gene ratio into a number i.e. calculate the ratio
#split the character values into two parts
split_result <- strsplit(GO_mrna$GeneRatio, "/")
split_result2 <- strsplit(GO_mrna_top$GeneRatio, "/")

#Convert the parts to numeric values and perform division
GeneRatioCalc <- sapply(split_result, function(x) as.numeric(x[1]) / as.numeric(x[2]))
GeneRatioCalc2 <- sapply(split_result2, function(x) as.numeric(x[1]) / as.numeric(x[2]))
#add the numeric result to a new column
GO_mrna$GeneRatioCalc <- GeneRatioCalc
GO_mrna_top$GeneRatioCalc2 <- GeneRatioCalc2

GO_mrna$log10_padj <- -log10(GO_mrna$p.adjust)
GO_mrna_top$log10_padj <- -log10(GO_mrna_top$p.adjust)
write.csv(GO_mrna, "E:/GO_mrna.csv", row.names = TRUE)
write.csv(GO_mrna_top, "E:/GO_mrna_top.csv", row.names = TRUE)

GO_mrna_top_up <- filter(GO_mrna_top, direction == "Up")
GO_mrna_top_down <- filter(GO_mrna_top, direction == "Down")

GO_mrna_top_up <- GO_mrna_top_up %>%
  arrange(ONTOLOGY, GeneRatioCalc2)


GO_mrna_top_up$group <- factor(GO_mrna_top_up$group, levels = c("sex-cords", "adenoma", "SBT", "HGSOC"))
GO_mrna_top_down$group <- factor(GO_mrna_top_down$group, levels = c("sex-cords", "adenoma", "SBT", "HGSOC"))
GO_mrna_top$direction <- factor(GO_mrna_top$direction, levels = c("Up", "Down"))
GO_mrna_top$group <- factor(GO_mrna_top$group, levels = c("sex-cords", "adenoma", "SBT", "HGSOC"))

GO_mrna_top_up_graph <- ggplot(GO_mrna_top_up, aes(x = group, y = Description)) + 
  geom_point(aes(color = log10_padj, size = Count)) + 
  facet_wrap(~ ONTOLOGY) + 
  scale_size_continuous() + 
  scale_color_gradient(low = "green", high = "purple") + 
  labs(x = "Group", y = "GO Term", 
       title = "Gene Ontology Enrichment Analysis", 
       subtitle = "Significantly Differentially Upregulated Genes (Group/Control)",
       color = "-log10(p.adj)", 
       size = "Gene Count") + 
  theme_gray() +
  theme(
    axis.text.x = element_text(angle = 90, size = 10.0, vjust = 0.5),
    axis.text.y = element_text(size = 10.0, vjust = 0.5),
    axis.title.x = element_text(size = 13.0, vjust = -3.0),
    axis.title.y = element_text(size = 13.0, vjust = 3.0),
    text = element_text(size = 10.0),
    plot.title = element_text(vjust = +3.0, hjust = 0.5),
    plot.margin = margin(1,1,1,1, "cm")
    )
 


GO_mrna_top_up_graph

ggsave("E:/GO_mrna_top_up_graph.png", plot = GO_mrna_top_up_graph, width = 12, height = 10, dpi = 600)

GO_mrna_top_down_graph <- ggplot(GO_mrna_top_down, aes(x = group, y = Description)) + 
  geom_point(aes(color = log10_padj, size = Count)) + 
  facet_wrap(~ ONTOLOGY) + 
  scale_size_continuous() + 
  scale_color_gradient(low = "green", high = "purple") + 
  labs(x = "Group", y = "GO Term", 
       title = "Gene Ontology Enrichment Analysis", 
       subtitle = "Significantly Differentially Downregulated Genes (Group/Control)",
       color = "-log10(p.adj)", 
       size = "Gene Count") + 
  theme_gray() +
  theme(
    axis.text.x = element_text(angle = 90, size = 10.0, vjust = 0.5),
    axis.text.y = element_text(size = 10.0, vjust = 0.5),
    axis.title.x = element_text(size = 13.0, vjust = -3.0),
    axis.title.y = element_text(size = 13.0, vjust = 3.0),
    text = element_text(size = 10.0),
    plot.title = element_text(vjust = +3.0, hjust = 0.5),
    plot.margin = margin(1,1,1,1, "cm")
  )

GO_mrna_top_down_graph
ggsave("E:/GO_mrna_top_down_graph.png", plot = GO_mrna_top_down_graph, width = 12, height = 10, dpi = 1200)

