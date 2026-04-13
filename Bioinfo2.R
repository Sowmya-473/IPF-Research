# Load GSE99621
gse99 <- read.table(
  "/Users/sowmi/Desktop/IPF_Research/GSE99621_datamatrix.txt",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

# Inspect
dim(gse99)
head(gse99[, 1:5])
colnames(gse99)

# Check count type
cat("\nMin:", min(gse99), "\n")
cat("Max:", max(gse99), "\n")
cat("Median:", median(as.matrix(gse99)), "\n")
any(gse99 != round(gse99))


# Install GEOquery to get sample metadata
if (!require("GEOquery", quietly = TRUE))
  BiocManager::install("GEOquery")

library(GEOquery)

# Download metadata for GSE99621
gse_meta <- getGEO("GSE99621", GSEMatrix = FALSE)

# Extract sample information
sample_info <- lapply(GSMList(gse_meta), function(x) {
  c(gsm = Meta(x)$geo_accession,
    title = Meta(x)$title,
    source = Meta(x)$source_name_ch1,
    characteristics = paste(Meta(x)$characteristics_ch1, collapse = "; "))
})

# Convert to dataframe
sample_df <- do.call(rbind, lapply(sample_info, as.data.frame))
sample_df <- as.data.frame(do.call(rbind, sample_info))

print(sample_df)


# Get the sample ID mapping from GEO
sample_mapping <- lapply(GSMList(gse_meta), function(x) {
  c(gsm       = Meta(x)$geo_accession,
    title     = Meta(x)$title,
    sample_id = Meta(x)$description)
})

mapping_df <- as.data.frame(do.call(rbind, sample_mapping))
print(mapping_df)

# Build metadata table matching OtB codes to conditions
metadata99 <- data.frame(
  sample_id = c("OtB5154","OtB5155","OtB5156","OtB5157","OtB5256",
                "OtB5159","OtB5160","OtB5161","OtB5162","OtB5163",
                "OtB5164","OtB5165","OtB5170","OtB5172","OtB5259",
                "OtB5178","OtB5179","OtB5180","OtB5166","OtB5254",
                "OtB5169","OtB5174","OtB5175","OtB5258","OtB5177",
                "OtB5264"),
  title     = c("HC11","HC12","HC21","HC22","HC23",
                "HC31","HC32","HC33","IPF1n1","IPF1n2",
                "IPF1n3","IPF1n4","IPF2n1","IPF2n2","IPF2n3",
                "IPF3n1","IPF3n2","IPF3n3","IPF1s1","IPF1s2",
                "IPF1s3","IPF2s1","IPF2s2","IPF2s3","IPF2s4",
                "IPF3s1"),
  condition = c(rep("Control", 8), rep("IPF", 18))
)

# Make sure column order in gse99 matches metadata
gse99_ordered <- gse99[, metadata99$sample_id]

# Confirm order matches
all(colnames(gse99_ordered) == metadata99$sample_id)

# Set rownames
rownames(metadata99) <- metadata99$sample_id

# Check
print(metadata99[, c("title", "condition")])


# Step 1: Create DESeq2 object
dds99 <- DESeqDataSetFromMatrix(
  countData = gse99_ordered,
  colData   = metadata99,
  design    = ~ condition
)

# Step 2: Filter low count genes
dds99 <- dds99[rowSums(counts(dds99)) >= 10, ]
cat("Genes after filtering:", nrow(dds99), "\n")

# Step 3: Run DESeq2
dds99 <- DESeq(dds99)
cat("DESeq2 finished!\n")

# Step 4: Extract results (IPF vs Control)
res99 <- results(dds99,
                 contrast = c("condition", "IPF", "Control"),
                 alpha = 0.05)

# Step 5: Summary
summary(res99)

# Step 6: Convert to dataframe
res99_df <- as.data.frame(res99)
res99_df$gene_symbol <- rownames(res99_df)

# Step 7: Filter significant DEGs
sig_degs99 <- subset(res99_df, padj < 0.05 & abs(log2FoldChange) > 1)
sig_degs99 <- sig_degs99[order(sig_degs99$padj), ]

# Results
cat("Total significant DEGs:", nrow(sig_degs99), "\n")
cat("Upregulated in IPF:", sum(sig_degs99$log2FoldChange > 0), "\n")
cat("Downregulated in IPF:", sum(sig_degs99$log2FoldChange < 0), "\n")

# Top 20
head(sig_degs99[, c("gene_symbol", "log2FoldChange", "padj")], 20)


# Save results
write.csv(sig_degs99,
          "/Users/sowmi/Desktop/IPF_Research/GSE99621_DEGs.csv",
          row.names = TRUE)

write.csv(res99_df,
          "/Users/sowmi/Desktop/IPF_Research/GSE99621_all_results.csv",
          row.names = TRUE)

cat("Saved!\n")

# ---- VOLCANO PLOT ----
library(EnhancedVolcano)

EnhancedVolcano(res99_df,
                lab      = res99_df$gene_symbol,
                x        = "log2FoldChange",
                y        = "padj",
                title    = "GSE99621: IPF vs Healthy Control",
                subtitle = "Lung Tissue RNA-seq (8 HC vs 18 IPF)",
                pCutoff  = 0.05,
                FCcutoff = 1.0,
                pointSize = 2,
                labSize  = 3,
                col      = c("grey70", "grey70", "steelblue", "firebrick"),
                drawConnectors = TRUE,
                widthConnectors = 0.3
)

ggsave("/Users/sowmi/Desktop/IPF_Research/volcano_GSE99621.pdf",
       width = 12, height = 10)

# ---- PCA PLOT ----
library(ggplot2)

vsd99 <- vst(dds99, blind = FALSE)

pcaData99 <- plotPCA(vsd99,
                     intgroup = "condition",
                     returnData = TRUE)

percentVar99 <- round(100 * attr(pcaData99, "percentVar"))

ggplot(pcaData99, aes(PC1, PC2, color = condition, label = name)) +
  geom_point(size = 4) +
  scale_color_manual(values = c("IPF" = "firebrick",
                                "Control" = "steelblue")) +
  xlab(paste0("PC1: ", percentVar99[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar99[2], "% variance")) +
  ggtitle("PCA: IPF vs Control (GSE99621)") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, size = 14))

ggsave("/Users/sowmi/Desktop/IPF_Research/PCA_GSE99621.pdf",
       width = 8, height = 6)

cat("Figures saved!\n")