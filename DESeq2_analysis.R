setwd("C:/Users/judo4/Desktop/OMICS_Analysis/OMICS1_Analysis")


# load libraries
library(DESeq2)
#BiocManager::install("DESeq2")
library(tidyverse)
library(dplyr)
#library(airway)

# Reprogrammed samples
SRR3414629 <- read.table("SRR3414629.counts", header = FALSE)
SRR3414630 <- read.table("SRR3414630.counts", header = FALSE)
SRR3414631 <- read.table("SRR3414631.counts", header = FALSE)

head(SRR3414629)

# control samples
SRR3414635 <- read.table("SRR3414635.counts", header = FALSE)
SRR3414636 <- read.table("SRR3414636.counts", header = FALSE)
SRR3414637 <- read.table("SRR3414637.counts", header = FALSE)

summary(SRR3414629)


# Remove the last rows containing .counts summary
# Define the rows to be removed
rows_to_remove <- c(57180, 57181, 57182, 57183, 57184)

# Remove the specified rows from all replicates
SRR3414629 <- SRR3414629[-rows_to_remove, ]
SRR3414630 <- SRR3414630[-rows_to_remove, ]
SRR3414631 <- SRR3414631[-rows_to_remove, ]

SRR3414635 <- SRR3414635[-rows_to_remove, ]
SRR3414636 <- SRR3414636[-rows_to_remove, ]
SRR3414637 <- SRR3414637[-rows_to_remove, ]

# Combine reprogrammed samples
merged_data <- full_join(SRR3414629, SRR3414630, by = "V1") %>%
  full_join(SRR3414631, by = "V1") %>%
  full_join(SRR3414635, by = "V1") %>%
  full_join(SRR3414636, by = "V1") %>%
  full_join(SRR3414637, by = "V1")

# Set the column names
colnames(merged_data) <- c("Gene_ID", "SRR3414629", "SRR3414630", "SRR3414631", "SRR3414635", "SRR3414636", "SRR3414637")

# Set the Gene_ID column as row names
rownames(merged_data) <- merged_data$Gene_ID

# Remove the Gene_ID column
combined_data <- merged_data[, -1]

summary(combined_data)

# Create a colData
# collect column names combined_data is your dataset
sample_names <- colnames(combined_data)

# Define metadata for each sample
colData <- data.frame(
  SampleID = sample_names,
  Samples = c("SRR3414629", "SRR3414630", "SRR3414631", "SRR3414635", "SRR3414636", "SRR3414637"),
  Condition = c("Reprogrammed", "Reprogrammed", "Reprogrammed", "Non_Reprogrammed", "Non_Reprogrammed", "Non_Reprogrammed")
  # Add more metadata columns as needed, such as batch information, treatment groups, etc.
)

# Convert SampleID to rownames
rownames(colData) <- colData$SampleID
colData$SampleID <- NULL  # Remove SampleID column

# Check the colData
head(colData)



# verifying dataset
# making sure the row names in colData matches to column names in combined_data
all(colnames(combined_data) %in% rownames(colData))

# are they in the same order
all (colnames(combined_data) == rownames(colData))


# Step 2: construct a DESeqDataSet object ----------

dds <- DESeqDataSetFromMatrix(countData = combined_data,
                              colData = colData,
                              design = ~ Condition)


summary(dds)


# pre-filtering: removing rows with low gene counts
# keeping rows that have at least 10 reads total
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds



# This is very important otherwise DESeq function will choose alphabetically
# set the factor level
dds$Condition <- relevel(dds$Condition, ref = "Non_Reprogrammed")

# The function DESeq runsin  the following functions in order:

# dds <- estimateSizeFactors(dds)
# dds <- estimateDispersions(dds)
# dds <- nbinomWaldTest(dds)

# NOTE: collapse technical replicates but this experiment uses biological replicate

# Step 3: Run DESeq ----------------------
dds <- DESeq(dds)
res <- results(dds)

res


#Dispersion plot
# This plot shows the relationship between the variation and mean
plotDispEsts(dds)

# Plots to visualize our dds before normalization.
# I had to combined with after normalization to show the difference.
# So DESeq2 should be performed before this step.
library(gridExtra)  # For arranging plots

# Compute PCA coordinates
pca <- prcomp(t(counts(dds)))

# Create dataframe with PCA coordinates and metadata for raw and normalized counts
pca_data <- as.data.frame(pca$x)
pca_data$Sample <- colData(dds)$Samples
pca_data$Condition <- colData(dds)$Condition

normalized_counts <- counts(dds, normalized = TRUE)
pca_deseq_norm <- prcomp(t(normalized_counts))
pca_data_deseq_norm <- as.data.frame(pca_deseq_norm$x)
pca_data_deseq_norm$Samples <- colData(dds)$Samples
pca_data_deseq_norm$Condition <- colData(dds)$Condition

# Get minimum and maximum values for PC1 and PC2 from raw counts
pc1_range <- c(min(pca_data$PC1), max(pca_data$PC1))
pc2_range <- c(min(pca_data$PC2), max(pca_data$PC2))

# Plot PCA with raw counts (without label)
plot_raw <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, shape = Sample)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1: ", round(summary(pca)$importance[2, 1], 2), "% variance"),
       y = paste0("PC2: ", round(summary(pca)$importance[2, 2], 2), "% variance"))+
  theme_minimal()+
  theme(legend.title = element_blank(), legend.position = "none")


# Plot PCA with DESeq2 normalized counts (using raw count scale)
plot_deseq_norm <- ggplot(pca_data_deseq_norm, aes(x = PC1, y = PC2, color = Condition, shape = Samples)) +
  geom_point(size = 3) +
  labs(x = paste0("PC1: ", round(summary(pca_deseq_norm)$importance[2, 1], 2), "% variance"),
       y = paste0("PC2: ", round(summary(pca_deseq_norm)$importance[2, 2], 2), "% variance")) +
  theme_minimal() +
  scale_x_continuous(limits = pc1_range) +  # Set scale based on raw count data
  scale_y_continuous(limits = pc2_range)  # Set scale based on raw count data

# Arrange plots side-by-side
grid.arrange(plot_raw, plot_deseq_norm, ncol = 2)






# Select counts for sample SRR3414635 from raw counts
raw_counts_SRR3414635 <- as.matrix(combined_data[, "SRR3414635"])

# Calculate the total counts for sample SRR3414635
total_counts_SRR3414635 <- sum(raw_counts_SRR3414635)

# Calculate the RPKM normalized counts
rpkm_normalized_counts_SRR3414635 <- (raw_counts_SRR3414635 * 1e6) / (total_counts_SRR3414635)

# Take the log10 of the raw counts
log2_raw_counts_SRR3414635 <- log2(raw_counts_SRR3414635)

# Take the log10 of the RPKM normalized counts
log2_rpkm_normalized_counts_SRR3414635 <- log2(rpkm_normalized_counts_SRR3414635)

# Select counts for sample SRR3414635 from normalized counts
normalized_counts_SRR3414635 <- log2(as.matrix(normalized_counts[, "SRR3414635"]))

# General title
main_title <- "Distribution of Raw and Normalized Counts (Log2) of Sample SRR3414635"

# Create separate plots for raw counts, normalized counts, and RPM normalized counts
plot_raw_counts <- ggplot(data.frame(Counts = log2_raw_counts_SRR3414635), aes(x = Counts)) +
  geom_histogram(binwidth = 0.5, fill = "blue", alpha = 0.7) +
  labs(title = NULL, x = "Reads Counts", y = "Frequency") +
  theme_minimal()

plot_normalized_counts <- ggplot(data.frame(Counts = normalized_counts_SRR3414635), aes(x = Counts)) +
  geom_histogram(binwidth = 0.5, fill = "red", alpha = 0.7) +
  labs(title = NULL, x = "Reads Counts", y = "Frequency") +
  theme_minimal()

plot_rpkm_normalized_counts <- ggplot(data.frame(Counts = log2_rpkm_normalized_counts_SRR3414635), aes(x = Counts)) +
  geom_histogram(binwidth = 0.5, fill = "green", alpha = 0.7) +
  labs(title = NULL, x = "Reads Counts", y = "Frequency") +
  theme_minimal()

# Print the plots side by side with a shared title and legend
grid.arrange(
  plot_raw_counts + ggtitle(main_title) + theme(legend.position = "none"),
  plot_normalized_counts + theme(legend.position = "none"),
  plot_rpkm_normalized_counts + theme(legend.position = "none"),
  ncol = 3
)

# Add legend separately
# Legend
Legend_Labels <- c("Raw Counts", "DESeq2 Normalized", "RPM Normalized")
legend_colors <- c("blue", "red", "green")

legend <- ggplot() +
  geom_point(aes(x = 1, y = 1, color = Legend_Labels), size = 3) +
  scale_color_manual(values = legend_colors, labels = Legend_Labels) +
  theme_void() +
  theme(legend.position = "bottom")

# Print the legend
legend














# I tried to use Lg2 Transformation to explain the distribution of padj
# Log 2 transformation
rld <- rlogTransformation(dds, blind = FALSE)

hist(assay(rld))

PCAA <- plotPCA(rld, intgroup="Condition")
PCAA <- geom_text(aes(label=name), size=2.5)
PCAA
# Histogram of p values for all tests. The area shaded in blue indicates the subset of those that pass the filtering, the area in khaki those that do not pass:
use <- res$baseMean > metadata(res)$filterThreshold
h1 <- hist(res$pvalue, breaks=0:50/50, plot=FALSE)
h2 <- hist(res$pvalue[use], breaks=0:50/50, plot=FALSE)
colori <- c(`do not pass`="khaki", `pass`="powderblue")

# Create the bar plot
barplot(height = rbind(h1$counts, h2$counts), beside = FALSE,
        col = colori, space = 0, main = "Histogram of padj for all tests",
        ylab = "Frequency")

# Add text labels for p-value ranges
text(x = c(0, length(h1$counts)), y = 0, label = paste(c(0, 1)),
     adj = c(0.5, 1.7), xpd = NA)

# Add legend
legend("topright", fill = rev(colori), legend = rev(names(colori)))



# Explore Results ----------------

summary(res)
res

# genes with padj < 0.01
res0.01 <- results(dds, alpha = 0.01)

summary(res0.01)

# sum of the number of genes with padj < 0.01
sum(res0.01$padj < 0.01, na.rm=TRUE)

res0.01

library("apeglm")
resLFC <- lfcShrink(dds, coef=2, res=res0.01)
resLFC


plotMA(resLFC)
plotMA(resLFC, ylim = c(-6.5,6.5))
topGene <- rownames(resLFC)[which.min(res$padj)]
with(res[topGene, ], {
  points(baseMean, log2FoldChange, col="red", cex=2, lwd=2)
  text(baseMean, log2FoldChange, topGene, pos=2, col="black")
})


# visualize the top Gene
topGene <- rownames(res)[which.min(res$padj)]
plotCounts(dds, gene = topGene, intgroup=c("Condition"))


# Normalized counts for a single gene over treatment group.
# 
# We can also make custom plots using the ggplot function from the ggplot2 package (figures below).
# install.packages("ggbeeswarm")
library("ggbeeswarm")
geneCounts <- plotCounts(dds, gene = topGene, intgroup = c("Condition","Samples"),
                         returnData = TRUE)
ggplot(geneCounts, aes(x = Condition, y = count, color = Condition, shape = Samples)) +
  scale_y_log10() +  geom_beeswarm(cex = 3)

# # Annotation
# install.packages("AnnotationDbi")
# install.packages("org.Hs.eg.db")
# 
# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Mm.eg.db")

library("AnnotationDbi")
library("org.Mm.eg.db")

columns(org.Mm.eg.db)

ens.str <- substr(rownames(res), 1, 18)
res$symbol <- mapIds(org.Mm.eg.db,
                     keys=ens.str,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
res$entrez <- mapIds(org.Mm.eg.db,
                     keys=ens.str,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")

resOrdered <- res[order(res$pvalue),]
head(resOrdered)

sum(resOrdered$padj < 0.01, na.rm=TRUE)


# Just incase the isoforms does not show try this
#keys(org.Mm.eg.db, keytype="ENSEMBLTRANS")


# checking and replacing all NAs with 1 for easy visualization of data
res$padj <- ifelse(is.na(res$padj), 1, res$padj)
# order differentially expressed genes by padj
resOrder <- res [order(res$padj), ]

# Define top differentially expressed genes
# Define threshold (adjust as needed)
pval_threshold <- 0.01  # Adjust threshold for significance (e.g., 0.01)

# Filter genes with adjusted p-value below threshold
top_DE_genes <- resOrder[resOrder$padj < pval_threshold, ]

# Plot a barchart for top 10 Genes
# Extract top 10 DE genes
top_10_DE_genes <- head(top_DE_genes, 10)

bar_width <- 0.89

# Create a bar plot for visualization
barplot(top_10_DE_genes$log2FoldChange, names.arg = top_10_DE_genes$symbol, las = 2,
        col = ifelse(top_10_DE_genes$log2FoldChange > 0, "blue", "red"),
        main = "Top 10 Differentially Expressed Genes",
        xlab = "Genes", ylab = "log2FoldChange", width = bar_width)

# Add adjusted p-values as text directly above each bar in exponential format
text(x = 1:10, y = top_10_DE_genes$log2FoldChange + 0.2, 
     labels = sprintf("%.1e", top_10_DE_genes$padj),
     pos = 3, cex = 0.8, xpd = TRUE)

# Increase width of each bin to accommodate larger p-values
par(mar=c(5, 8, 4, 2))  # Adjust the margin to make space for longer x-axis labels


# Add a legend for log2FoldChange direction
legend("topright", legend = c("Upregulated", "Downregulated"),
       fill = c("blue", "red"), bg = "white", inset = 0.05)



# Create data.frame with desired information
OMICS1_DE_genes <- data.frame(Gene_ID = rownames(top_DE_genes),
                              log2FoldChange = top_DE_genes$log2FoldChange,
                              Gene_Annotation = top_DE_genes$symbol,
                              padj = top_DE_genes$padj)

# Write data.frame to CSV file
write.csv(OMICS1_DE_genes, file = "DE_OMICS1.csv", row.names = FALSE)

print("Top differentially expressed genes written to DE_OMICS1.csv")

summary(top_DE_genes)



# PATHWAY ENRICHMENT ANALYSIS
# Create data.frame with desired information for enrichment analysis
OMICS1_DE_res <- data.frame(Gene_Annotation = res$symbol,
                            log2FoldChange = res$log2FoldChange,
                            pvalue = res$pvalue,
                            padj = res$padj)


# Write dataframe to CSV file for enrichment analysis
write.csv(OMICS1_DE_res, file = "OMICS1_DE_res.csv", row.names = FALSE)