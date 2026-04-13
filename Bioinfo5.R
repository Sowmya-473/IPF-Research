# Load GSE310473 processed data
gse473 <- read.table(
  "/Users/sowmi/Desktop/IPF_Research/GSE310473_processed.txt.gz",
  header = TRUE,
  sep = "\t",
  row.names = 1
)

# Inspect
dim(gse473)
head(gse473[, 1:6])
colnames(gse473)


# ---- SEPARATE COUNT MATRICES ----

# ION363 comparison (ASO drug vs scrambled)
counts_ion <- gse473[, c("scr.ION.363_1", "scr.ION.363_2", 
                         "scr.ION.363_3", "scr.ION.363_4",
                         "ION.363_1",     "ION.363_2", 
                         "ION.363_3",     "ION.363_4")]

# siFUS comparison (siRNA knockdown vs scrambled)
counts_sifus <- gse473[, c("scr.siFUS_1", "scr.siFUS_2",
                           "scr.siFUS_3", "scr.siFUS_4",
                           "siFUS_1",    "siFUS_2",
                           "siFUS_3",    "siFUS_4")]

# Fix NA values
counts_ion[is.na(counts_ion)]     <- 0
counts_sifus[is.na(counts_sifus)] <- 0

# Check for decimals
cat("ION363 decimals:", any(counts_ion != round(counts_ion)), "\n")
cat("siFUS decimals:", any(counts_sifus != round(counts_sifus)), "\n")

# Gene info
gene_info473 <- gse473[, c("Gene.Symbol", "Gene.Biotype", "Gene.Description")]

# ---- METADATA ----
meta_ion <- data.frame(
  sample    = colnames(counts_ion),
  condition = c(rep("Scrambled", 4), rep("ION363", 4))
)
rownames(meta_ion) <- meta_ion$sample

meta_sifus <- data.frame(
  sample    = colnames(counts_sifus),
  condition = c(rep("Scrambled", 4), rep("siFUS", 4))
)
rownames(meta_sifus) <- meta_sifus$sample

cat("\nION363 metadata:\n")
print(meta_ion[, "condition", drop=FALSE])
cat("\nsiFUS metadata:\n")
print(meta_sifus[, "condition", drop=FALSE])


library(DESeq2)

# ---- DESeq2: ION363 vs Scrambled ----
dds_ion <- DESeqDataSetFromMatrix(
  countData = counts_ion,
  colData   = meta_ion,
  design    = ~ condition
)
dds_ion <- dds_ion[rowSums(counts(dds_ion)) >= 10, ]
dds_ion <- DESeq(dds_ion)

res_ion <- results(dds_ion,
                   contrast = c("condition", "ION363", "Scrambled"),
                   alpha = 0.05)

# ---- DESeq2: siFUS vs Scrambled ----
dds_sifus <- DESeqDataSetFromMatrix(
  countData = counts_sifus,
  colData   = meta_sifus,
  design    = ~ condition
)
dds_sifus <- dds_sifus[rowSums(counts(dds_sifus)) >= 10, ]
dds_sifus <- DESeq(dds_sifus)

res_sifus <- results(dds_sifus,
                     contrast = c("condition", "siFUS", "Scrambled"),
                     alpha = 0.05)

# ---- Summaries ----
cat("=== ION363 vs Scrambled ===\n")
summary(res_ion)

cat("\n=== siFUS vs Scrambled ===\n")
summary(res_sifus)


# ---- Extract DEGs from both ----
res_ion_df <- as.data.frame(res_ion)
res_ion_df$ensembl_id  <- rownames(res_ion_df)
res_ion_df$gene_symbol <- gene_info473[rownames(res_ion_df), "Gene.Symbol"]

res_sifus_df <- as.data.frame(res_sifus)
res_sifus_df$ensembl_id  <- rownames(res_sifus_df)
res_sifus_df$gene_symbol <- gene_info473[rownames(res_sifus_df), "Gene.Symbol"]

# Significant DEGs
sig_ion <- subset(res_ion_df, 
                  padj < 0.05 & abs(log2FoldChange) > 1)
sig_ion <- sig_ion[order(sig_ion$padj), ]

sig_sifus <- subset(res_sifus_df,
                    padj < 0.05 & abs(log2FoldChange) > 1)
sig_sifus <- sig_sifus[order(sig_sifus$padj), ]

cat("ION363 DEGs (|LFC|>1):", nrow(sig_ion), "\n")
cat("siFUS DEGs (|LFC|>1):", nrow(sig_sifus), "\n")

# Top 10 ION363 DEGs
cat("\nTop 10 ION363 DEGs:\n")
print(head(sig_ion[, c("gene_symbol", "log2FoldChange", "padj")], 10))

# Top 10 siFUS DEGs
cat("\nTop 10 siFUS DEGs:\n")
print(head(sig_sifus[, c("gene_symbol", "log2FoldChange", "padj")], 10))

# ---- THREE-WAY OVERLAP ----
# Genes changed by BOTH FUS silencing methods AND dysregulated in IPF disease
genes_ion   <- sig_ion$gene_symbol[!is.na(sig_ion$gene_symbol)]
genes_sifus <- sig_sifus$gene_symbol[!is.na(sig_sifus$gene_symbol)]

# Two-way: ION363 ∩ siFUS (genes changed by BOTH silencing methods)
overlap_ion_sifus <- intersect(genes_ion, genes_sifus)
cat("\nION363 ∩ siFUS overlap:", length(overlap_ion_sifus), "genes\n")

# Three-way: ION363 ∩ siFUS ∩ IPF disease (GSE99621)
overlap_3way <- intersect(overlap_ion_sifus, genes_99)
cat("Three-way overlap (ION363 ∩ siFUS ∩ IPF):", 
    length(overlap_3way), "genes\n")
print(sort(overlap_3way))

# Save everything
write.csv(sig_ion,
          "/Users/sowmi/Desktop/IPF_Research/GSE310473_ION363_DEGs.csv",
          row.names = TRUE)
write.csv(sig_sifus,
          "/Users/sowmi/Desktop/IPF_Research/GSE310473_siFUS_DEGs.csv",
          row.names = TRUE)
write.csv(data.frame(gene = overlap_3way),
          "/Users/sowmi/Desktop/IPF_Research/overlap_3way.csv",
          row.names = FALSE)

cat("\nAll saved!\n")



library(VennDiagram)

# ---- EXTENDED OVERLAP ANALYSIS ----

# Also check: how many of your original 122 genes appear in ION363 dataset?
overlap_pcls_ion <- intersect(overlap_genes, genes_ion)
cat("Original 122 ∩ ION363:", length(overlap_pcls_ion), "genes\n")

# How many of original 122 appear in siFUS?
overlap_pcls_sifus <- intersect(overlap_genes, genes_sifus)
cat("Original 122 ∩ siFUS:", length(overlap_pcls_sifus), "genes\n")

# Four-way: PCLS ∩ ION363 ∩ siFUS ∩ IPF disease
overlap_4way <- intersect(overlap_genes, 
                          intersect(overlap_ion_sifus, genes_99))
cat("Four-way overlap (PCLS ∩ ION363 ∩ siFUS ∩ IPF):", 
    length(overlap_4way), "genes\n")
print(sort(overlap_4way))

# Install scales if needed
if (!require("scales")) install.packages("scales")
library(scales)

# Now run the Venn diagram
venn.diagram(
  x = list(
    "ION363\n(GSE310473)"     = genes_ion,
    "siFUS\n(GSE310473)"      = genes_sifus,
    "IPF disease\n(GSE99621)" = genes_99
  ),
  filename   = "/Users/sowmi/Desktop/IPF_Research/venn_3way_473.png",
  output     = TRUE,
  imagetype  = "png",
  height     = 2000,
  width      = 2500,
  resolution = 300,
  col  = c("firebrick", "steelblue", "darkgreen"),
  fill = c(alpha("firebrick",  0.3),
           alpha("steelblue",  0.3),
           alpha("darkgreen",  0.3)),
  cex          = 1.5,
  fontface     = "bold",
  cat.cex      = 1.1,
  cat.fontface = "bold",
  margin       = 0.12,
  main         = "FUS silencing methods vs IPF disease genes",
  main.cex     = 1.4,
  main.fontface = "bold"
)

cat("Venn diagram saved!\n")

# Also print the overlap numbers
cat("\nOriginal 122 ∩ ION363:", length(overlap_pcls_ion), "genes\n")
cat("Original 122 ∩ siFUS:", length(overlap_pcls_sifus), "genes\n")
cat("Four-way overlap:", length(overlap_4way), "genes\n")
print(sort(overlap_4way))

# ---- MULTI-DATASET SUMMARY BAR CHART ----
library(ggplot2)

summary_df <- data.frame(
  dataset     = c("GSE310475\n(PCLS ION363)",
                  "GSE310473\n(Fibroblast ION363)",
                  "GSE310473\n(Fibroblast siFUS)",
                  "GSE99621\n(IPF vs Healthy)",
                  "ION363 ∩ siFUS\n(both methods)",
                  "3-way overlap\n(both methods + IPF)",
                  "Original 122\n∩ ION363 fibroblast",
                  "Original 122\n∩ siFUS fibroblast"),
  n_genes     = c(633, 1670, 115, 1486, 87, 7,
                  length(overlap_pcls_ion),
                  length(overlap_pcls_sifus)),
  category    = c("Primary", "Validation", "Validation",
                  "Disease", "Key", "Key", "Key", "Key")
)

ggplot(summary_df, 
       aes(x = reorder(dataset, n_genes), 
           y = n_genes, 
           fill = category)) +
  geom_bar(stat = "identity", width = 0.7) +
  geom_text(aes(label = n_genes), 
            hjust = -0.2, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = c(
    "Primary"    = "steelblue",
    "Validation" = "coral",
    "Disease"    = "darkgreen",
    "Key"        = "firebrick"
  )) +
  coord_flip() +
  expand_limits(y = max(summary_df$n_genes) * 1.15) +
  labs(
    title    = "Multi-dataset FUS-IPF transcriptomic analysis",
    subtitle = "Number of significant DEGs per comparison",
    x = NULL,
    y = "Number of DEGs",
    fill = "Category"
  ) +
  theme_classic() +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 13),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.y   = element_text(size = 9)
  )

ggsave("/Users/sowmi/Desktop/IPF_Research/multidataset_summary.pdf",
       width = 10, height = 7)

cat("Saved!\n")

# Also print the overlap numbers we need
cat("\nOriginal 122 ∩ ION363 fibroblasts:", length(overlap_pcls_ion), "\n")
cat("Original 122 ∩ siFUS fibroblasts:", length(overlap_pcls_sifus), "\n")
cat("Four-way overlap:", length(overlap_4way), "\n")
cat("Four-way genes:", paste(sort(overlap_4way), collapse=", "), "\n")