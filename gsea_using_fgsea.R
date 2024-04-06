# Setting up environment ===================================================
setwd("C:/Users/judo4/Desktop/OMICS_Analysis/OMICS1_Analysis")

# Clean environment
rm(list = ls(all.names = TRUE)) # will clear all objects including hidden objects
gc() # free up memory and report the memory usage
options(max.print = .Machine$integer.max, scipen = 999, stringsAsFactors = F, dplyr.summarise.inform = F) # avoid truncated output in R console and scientific notation


# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(pheatmap)
library(clusterProfiler) # for PEA analysis
library('org.Mm.eg.db')
library(DOSE)
library(enrichplot) # for visualisations
library(ggupset) # for visualisations
library(fgsea)

# Set input path
in_path <- "data/" # input path, where your data is located
out_path <- "results/" # output path, where you want your results exported to
bg_path <- "data/" # folder with your background genes used for PEA



# # For data management
# install.packages('tidyverse')
# BiocManager::install("clusterProfiler")
# BiocManager::install("org.Hs.eg.db")
# # For visualisation
# install.packages('pheatmap')
# install.packages("DOSE")
# install.packages("enrichplot")
# install.packages("ggupset")
# BiocManager::install("fgsea")


# GENE SET ENRICHMENT ANALYSIS

# Steps for fgseaMultilel function

# Function: Adjacency matrix to list -------------------------
matrix_to_list <- function(pws){
  pws.l <- list()
  for (pw in colnames(pws)) {
    pws.l[[pw]] <- rownames(pws)[as.logical(pws[, pw])]
  }
  return(pws.l)
}

# function prepare_gmt
prepare_gmt <- function(gmt_file, genes_in_data, savefile = FALSE) {
  # Read in gmt file
  gmt <- gmtPathways(gmt_file)
  
  # Get unique genes from the gmt file
  hidden <- unique(unlist(gmt))
  
  # Convert gmt file to a matrix with genes as rows and annotations as columns
  mat <- matrix(NA, dimnames = list(hidden, names(gmt)),
                nrow = length(hidden), ncol = length(gmt))
  
  for (i in 1:dim(mat)[2]) {
    mat[, i] <- as.numeric(hidden %in% gmt[[i]])
  }
  
  # Subset to genes present in our data to avoid bias
  hidden1 <- intersect(genes_in_data, hidden)
  mat <- mat[hidden1, colnames(mat)[which(colSums(mat[hidden1, ]) > 5)]]
  
  # Convert matrix back to list
  final_list <- matrix_to_list(mat)
  
  # Save file if specified
  if (savefile) {
    saveRDS(final_list, file = paste0(gsub('.gmt', '', gmt_file), '_subset_', format(Sys.time(), '%d%m'), '.RData'))
  }
  
  print('Wohoo! .gmt conversion successful! :)')
  return(final_list)
}

# Analysis ====================================================
## 1. Read in data -----------------------------------------------------------
list.files(in_path)
df <- read.csv(paste0(in_path, 'OMICS1_DE_res.csv'), row.names = NULL)
df <- na.omit(df)
## 2. Prepare background genes -----------------------------------------------
# Download gene sets .gmt files
#https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp
# For GSEA
# Filter out the gmt files for KEGG, Reactome and GOBP
genes_in_data <- df$Gene_Annotation
list.files(bg_path)
gmt_files <- "data/mh.all.v2023.2.Mm.symbols.gmt"
gmt_files
bg_genes <- prepare_gmt(gmt_files, genes_in_data, savefile = FALSE)

# Ranking using sign(df$log2fc)*(-log10(df$pval))) as a ranking metric
rankings <- sign(df$log2FoldChange)*(-log10(df$padj)) # we will use the signed p values from spatial DGE as ranking
names(rankings) <- df$Gene_Annotation # genes as names#
head(rankings)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)
max(rankings)
min(rankings)

# Some genes have such low p values that the signed pval is +- inf, we need to change it to the maximum * constant to avoid problems with fgsea
max_ranking <- max(rankings[is.finite(rankings)])
min_ranking <- min(rankings[is.finite(rankings)])
rankings <- replace(rankings, rankings > max_ranking, max_ranking * 1)
rankings <- replace(rankings, rankings < min_ranking, min_ranking * 1)
rankings <- sort(rankings, decreasing = TRUE) # sort genes by ranking
plot(rankings)


ggplot(data.frame(gene_symbol = names(rankings)[1:50], ranks = rankings[1:50]), aes(gene_symbol, ranks)) + 
  geom_point() +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))


## 4. Run GSEA ---------------------------------------------------------------
# Easy peasy! Run fgsea with the pathways 
GSEAres <- fgsea(pathways = bg_genes, # List of gene sets to check
                 stats = rankings,
                 scoreType = 'std', # in this case we have both pos and neg rankings. if only pos or neg, set to 'pos', 'neg'
                 #minSize = 10,
                 #maxSize = 500,
                 nproc = 1) # for parallelisation


head(GSEAres[order(padj)])

# Number of significant pathways at padj < 0.01
sum(GSEAres[, padj < 0.05])
sum(GSEAres[, pval < 0.05])

# Sorting the data by NES
GSEAres <- GSEAres[order(GSEAres$NES), ]

# Creating the bar chart
ggplot(GSEAres, aes(x = reorder(pathway, NES), y = NES, fill = padj < 0.05)) +
  geom_col() +  # Use geom_col() instead of geom_bar() to avoid 'stat' error
  scale_fill_manual(name = "Significance", values = c("TRUE" = "orange", "FALSE" = "gray"),
                    labels = c("Not Significant", "Significant")) +  # Assign colors and labels
  coord_flip() +  # Flip the coordinates to make it vertical
  theme_minimal() +  # Choose a minimal theme for the plot
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +  # Rotate x-axis labels
  labs(x = "Pathway", y = "NES") +  # Adding labels to x and y axes
  guides(fill = guide_legend(title = "padj"))  # Add legend title