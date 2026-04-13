install.packages(c("tidyverse", "pheatmap", "ggplot2"))
install.packages("BiocManager")
install.packages("EnhancedVolcano")
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "STRINGdb"))
BiocManager::install("EnhancedVolcano")
if (!require("BiocManager")) install.packages("BiocManager")

BiocManager::install(version = "3.22")
BiocManager::install(c("S4Vectors", "DESeq2"), force = TRUE)
# ---- LOOK AT GSE310475 FIRST ----

# Read the processed file
gse <- read.table(
  "/Users/sowmi/Desktop/IPF_Research/GSE310475_processed.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1
)
gse
# Basic inspection
dim(gse)        # How many genes x how many samples?
head(gse)       # Show first 6 rows
colnames(gse)   # What are the sample names?

# ---- CLEAN UP THE DATA ----

# Keep only the count columns (the 6 samples)
counts_gse <- gse[, c("Scr_1", "Scr_2", "Scr_3", "ASO_1", "ASO_2", "ASO_3")]

# Keep gene symbols as a reference
gene_info <- gse[, c("Gene.Symbol", "Gene.Biotype", "Gene.Description")]

# Check it looks right
head(counts_gse)
dim(counts_gse)   # should say 15803 x 6

# Check for any decimal numbers (important!)
# If TRUE, numbers are decimals. If FALSE, they are whole numbers.
any(counts_gse != round(counts_gse))

# How many missing values are there?
sum(is.na(counts_gse))

# Which rows have missing values? (show first few)
head(which(rowSums(is.na(counts_gse)) > 0))

# What does one of those rows look like?
counts_gse[which(rowSums(is.na(counts_gse)) > 0)[1], ]

# Replace all NA with 0
counts_gse[is.na(counts_gse)] <- 0

# Confirm it's fixed — should now say FALSE
any(counts_gse != round(counts_gse))

# Confirm no more NAs — should say 0
sum(is.na(counts_gse))

# Load DESeq2
library(DESeq2)

# Step 1: Create metadata table (tells DESeq2 which samples are which group)
metadata <- data.frame(
  sample    = c("Scr_1", "Scr_2", "Scr_3", "ASO_1", "ASO_2", "ASO_3"),
  condition = c("Scrambled", "Scrambled", "Scrambled", "ASO", "ASO", "ASO")
)
rownames(metadata) <- metadata$sample

# Check it looks right
print(metadata)

# Step 2: Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts_gse,
  colData   = metadata,
  design    = ~ condition
)

# Step 3: Filter out genes with very low counts
# (genes with total count < 10 across all samples are unreliable)
dds <- dds[rowSums(counts(dds)) >= 10, ]

# How many genes remain after filtering?
cat("Genes after filtering:", nrow(dds), "\n")

# Step 4: Run DESeq2
dds <- DESeq(dds)

cat("DESeq2 finished successfully!\n")

# Step 5: Extract results (ASO vs Scrambled)
res <- results(dds, 
               contrast = c("condition", "ASO", "Scrambled"),
               alpha = 0.05)

# Step 6: Summary of results
summary(res)

# Step 7: Convert to dataframe and add gene symbols
res_df <- as.data.frame(res)
res_df$ensembl_id  <- rownames(res_df)
res_df$gene_symbol <- gene_info[rownames(res_df), "Gene.Symbol"]

# Step 8: Filter significant DEGs
# Criteria: adjusted p-value < 0.05 AND absolute fold change > 1 (i.e. 2x difference)
sig_degs <- subset(res_df, padj < 0.05 & abs(log2FoldChange) > 1)

# Sort by adjusted p-value (most significant first)
sig_degs <- sig_degs[order(sig_degs$padj), ]

# How many significant DEGs?
cat("Total significant DEGs:", nrow(sig_degs), "\n")
cat("Upregulated (higher in ASO):", sum(sig_degs$log2FoldChange > 0), "\n")
cat("Downregulated (lower in ASO):", sum(sig_degs$log2FoldChange < 0), "\n")

# Show top 10 most significant genes
head(sig_degs[, c("gene_symbol", "log2FoldChange", "padj")], 10)

# ---- SAVE RESULTS ----

# Save all significant DEGs to CSV
write.csv(sig_degs, 
          "/Users/sowmi/Desktop/IPF_Research/GSE310475_DEGs.csv",
          row.names = TRUE)

# Save full results (all genes, not just significant)
write.csv(res_df,
          "/Users/sowmi/Desktop/IPF_Research/GSE310475_all_results.csv",
          row.names = TRUE)

cat("Files saved!\n")

# ---- MAKE VOLCANO PLOT ----
library(EnhancedVolcano)

EnhancedVolcano(res_df,
                lab     = res_df$gene_symbol,
                x       = "log2FoldChange",
                y       = "padj",
                title   = "GSE310475: ASO (FUS-silenced) vs Scrambled",
                subtitle = "IPF Precision-Cut Lung Slices",
                pCutoff = 0.05,
                FCcutoff = 1.0,
                pointSize = 2,
                labSize = 3,
                col = c("grey70", "grey70", "steelblue", "firebrick"),
                legendLabels = c("NS", "LFC only", "p-value only", "Significant"),
                drawConnectors = TRUE,
                widthConnectors = 0.3
)

# Save the plot
ggsave("/Users/sowmi/Desktop/IPF_Research/volcano_GSE310475.pdf",
       width = 12, height = 10)

cat("Volcano plot saved!\n")

library(pheatmap)

# Get variance stabilized data for visualization
vsd <- vst(dds, blind = FALSE)

# Get top 50 most significant DEGs
top50 <- head(sig_degs, 50)

# Extract their expression values
heat_mat <- assay(vsd)[rownames(top50), ]

# Replace Ensembl IDs with gene symbols as row names
rownames(heat_mat) <- top50$gene_symbol

# Create annotation bar showing which samples are which group
anno_col <- data.frame(
  Group = c("Scrambled", "Scrambled", "Scrambled", "ASO", "ASO", "ASO")
)
rownames(anno_col) <- colnames(heat_mat)

# Define colors
anno_colors <- list(
  Group = c(Scrambled = "steelblue", ASO = "firebrick")
)

# Plot heatmap
pheatmap(heat_mat,
         annotation_col  = anno_col,
         annotation_colors = anno_colors,
         scale           = "row",
         show_rownames   = TRUE,
         show_colnames   = TRUE,
         fontsize_row    = 7,
         fontsize_col    = 10,
         main            = "Top 50 DEGs: ASO vs Scrambled (GSE310475)",
         color           = colorRampPalette(c("steelblue", "white", "firebrick"))(100)
)

cat("Heatmap saved!\n")

library(ggplot2)

# Get PCA data
pcaData <- plotPCA(vsd, 
                   intgroup = "condition", 
                   returnData = TRUE)

# Get variance explained
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Plot
ggplot(pcaData, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 5) +
  geom_text(vjust = -1, size = 3.5) +
  scale_color_manual(values = c("ASO" = "firebrick", 
                                "Scrambled" = "steelblue")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: ASO vs Scrambled (GSE310475)") +
  theme_classic() +
  theme(legend.title = element_text(size = 12),
        legend.text  = element_text(size = 11),
        plot.title   = element_text(hjust = 0.5, size = 14))

ggsave("/Users/sowmi/Desktop/IPF_Research/PCA_GSE310475.pdf",
       width = 8, height = 6)

cat("PCA plot saved!\n")

