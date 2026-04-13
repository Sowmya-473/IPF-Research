install.packages("VennDiagram")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
# ---- OVERLAP ANALYSIS ----

# Get gene names from each dataset
genes_gse <- sig_degs$gene_symbol        # FUS-regulated genes (GSE310475)
genes_99  <- sig_degs99$gene_symbol      # IPF disease genes (GSE99621)

# Remove NAs
genes_gse <- genes_gse[!is.na(genes_gse)]
genes_99  <- genes_99[!is.na(genes_99)]

cat("GSE310475 DEGs:", length(genes_gse), "\n")
cat("GSE99621 DEGs:", length(genes_99), "\n")

# Find overlap
overlap_genes <- intersect(genes_gse, genes_99)
cat("Overlapping genes:", length(overlap_genes), "\n")
print(sort(overlap_genes))

library(VennDiagram)
library(ggplot2)

# Save overlap gene list
overlap_df <- data.frame(
  gene_symbol = overlap_genes,
  in_GSE310475 = overlap_genes %in% genes_gse,
  in_GSE99621  = overlap_genes %in% genes_99
)

# Add fold change info from both datasets
overlap_df$LFC_ASO_vs_Scr <- sig_degs$log2FoldChange[
  match(overlap_genes, sig_degs$gene_symbol)]

overlap_df$LFC_IPF_vs_Ctrl <- sig_degs99$log2FoldChange[
  match(overlap_genes, sig_degs99$gene_symbol)]

overlap_df$padj_GSE310475 <- sig_degs$padj[
  match(overlap_genes, sig_degs$gene_symbol)]

overlap_df$padj_GSE99621 <- sig_degs99$padj[
  match(overlap_genes, sig_degs99$gene_symbol)]

# Sort by significance
overlap_df <- overlap_df[order(overlap_df$padj_GSE310475), ]

# Save
write.csv(overlap_df,
          "/Users/sowmi/Desktop/IPF_Research/overlap_122_genes.csv",
          row.names = FALSE)

cat("Overlap saved!\n")

# ---- VENN DIAGRAM ----
venn.diagram(
  x = list(
    "IPF Disease\n(GSE99621)"        = genes_99,
    "FUS-regulated\n(GSE310475)"     = genes_gse
  ),
  filename = "/Users/sowmi/Desktop/IPF_Research/venn_diagram_final.png",
  output = TRUE,
  imagetype = "png",
  height = 2000,
  width  = 3000,
  resolution = 300,
  margin = 0.15,
  col  = c("steelblue", "firebrick"),
  fill = c(alpha("steelblue", 0.3), alpha("firebrick", 0.3)),
  cex = 1.8,
  fontface = "bold",
  cat.cex = 1.1,
  cat.fontface = "bold",
  cat.col = c("steelblue", "firebrick"),
  cat.pos = c(-20, 20),
  cat.dist = c(0.05, 0.05),
  main = "FUS-regulated vs IPF Disease Genes",
  main.cex = 1.5,
  main.fontface = "bold"
)

library(clusterProfiler)
library(org.Hs.eg.db)

# Convert gene symbols to Entrez IDs
overlap_entrez <- bitr(overlap_genes,
                       fromType = "SYMBOL",
                       toType   = "ENTREZID",
                       OrgDb    = org.Hs.eg.db)

cat("Genes successfully converted:", nrow(overlap_entrez), "\n")

# GO Biological Process enrichment
go_results <- enrichGO(
  gene          = overlap_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

# KEGG pathway enrichment
kegg_results <- enrichKEGG(
  gene          = overlap_entrez$ENTREZID,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

cat("GO terms found:", nrow(go_results), "\n")
cat("KEGG pathways found:", nrow(kegg_results), "\n")

# Plot GO results
barplot(go_results, 
        showCategory = 15,
        title = "GO Biological Processes — 122 Overlap Genes")

# Plot KEGG results  
barplot(kegg_results,
        showCategory = 15,
        title = "KEGG Pathways — 122 Overlap Genes")
#-----------------------------------------------------------------------------
# Fix 1: Use dotplot instead of barplot (works better with clusterProfiler v4)
dotplot(go_results, showCategory = 15,
        title = "GO Biological Processes — 122 Overlap Genes")

dotplot(kegg_results, showCategory = 15,
        title = "KEGG Pathways — 122 Overlap Genes")

# Fix 2: Also run enrichment on ALL GSE310475 DEGs (633 genes)
# This gives richer pathway results for your paper
all_gse_entrez <- bitr(genes_gse,
                       fromType = "SYMBOL",
                       toType   = "ENTREZID",
                       OrgDb    = org.Hs.eg.db)

go_gse <- enrichGO(
  gene          = all_gse_entrez$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  readable      = TRUE
)

kegg_gse <- enrichKEGG(
  gene          = all_gse_entrez$ENTREZID,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

cat("GO terms (all GSE310475 DEGs):", nrow(go_gse), "\n")
cat("KEGG pathways (all GSE310475 DEGs):", nrow(kegg_gse), "\n")

dotplot(go_gse, showCategory = 20,
        title = "GO Biological Processes — GSE310475 DEGs")

dotplot(kegg_gse, showCategory = 20,
        title = "KEGG Pathways — GSE310475 DEGs")

# Save all 4 enrichment plots as PDFs

# 1. GO overlap genes
pdf("/Users/sowmi/Desktop/IPF_Research/GO_overlap_122.pdf", width=10, height=8)
dotplot(go_results, showCategory=15,
        title="GO Biological Processes — 122 Overlap Genes")
dev.off()

# 2. KEGG overlap genes
pdf("/Users/sowmi/Desktop/IPF_Research/KEGG_overlap_122.pdf", width=10, height=8)
dotplot(kegg_results, showCategory=15,
        title="KEGG Pathways — 122 Overlap Genes")
dev.off()

# 3. GO all GSE310475 DEGs
pdf("/Users/sowmi/Desktop/IPF_Research/GO_GSE310475.pdf", width=12, height=10)
dotplot(go_gse, showCategory=20,
        title="GO Biological Processes — GSE310475 DEGs")
dev.off()

# 4. KEGG all GSE310475 DEGs
pdf("/Users/sowmi/Desktop/IPF_Research/KEGG_GSE310475.pdf", width=10, height=8)
dotplot(kegg_gse, showCategory=20,
        title="KEGG Pathways — GSE310475 DEGs")
dev.off()

cat("All plots saved!\n")

# Also run enrichment on GSE99621 DEGs for comparison
all_99_entrez <- bitr(genes_99,
                       fromType = "SYMBOL",
                       toType   = "ENTREZID",
                       OrgDb    = org.Hs.eg.db)

kegg_99 <- enrichKEGG(
  gene          = all_99_entrez$ENTREZID,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05
)

cat("KEGG pathways GSE99621:", nrow(kegg_99), "\n")

pdf("/Users/sowmi/Desktop/IPF_Research/KEGG_GSE99621.pdf", width=10, height=8)
dotplot(kegg_99, showCategory=20,
        title="KEGG Pathways — GSE99621 IPF vs Control")
dev.off()

# Save all enrichment results as CSV
write.csv(as.data.frame(go_gse),
          "/Users/sowmi/Desktop/IPF_Research/GO_GSE310475_results.csv")
write.csv(as.data.frame(kegg_gse),
          "/Users/sowmi/Desktop/IPF_Research/KEGG_GSE310475_results.csv")
write.csv(as.data.frame(kegg_results),
          "/Users/sowmi/Desktop/IPF_Research/KEGG_overlap_results.csv")

cat("All enrichment results saved!\n")

# Load the STRING interaction file
string_data <- read.table(
  "/Users/sowmi/Desktop/IPF_Research/string_interactions.tsv",
  header = TRUE,
  sep = "\t"
)

head(string_data)
colnames(string_data)

# Count connections for each gene (degree = number of interactions)
all_nodes <- c(string_data[,1], string_data[,2])
degree_table <- sort(table(all_nodes), decreasing = TRUE)

# Show top 20 hub genes
cat("Top 20 Hub Genes by Number of Connections:\n")
print(head(degree_table, 20))


library(ggplot2)

# Create hub gene dataframe
hub_df <- data.frame(
  gene   = names(head(degree_table, 15)),
  degree = as.numeric(head(degree_table, 15))
)

# Sort for plotting
hub_df <- hub_df[order(hub_df$degree, decreasing = TRUE), ]
hub_df$gene <- factor(hub_df$gene, levels = rev(hub_df$gene))

# Plot
ggplot(hub_df, aes(x = gene, y = degree, fill = degree)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "steelblue", high = "firebrick") +
  coord_flip() +
  labs(
    title = "Top 15 Hub Genes — PPI Network",
    subtitle = "122 FUS-regulated IPF Overlap Genes (STRING)",
    x = "Gene",
    y = "Number of Interactions (Degree)"
  ) +
  theme_classic() +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.y   = element_text(size = 11, face = "bold")
  ) +
  guides(fill = "none")

ggsave("/Users/sowmi/Desktop/IPF_Research/hub_genes_barplot.pdf",
       width = 8, height = 7)

cat("Hub gene plot saved!\n")

# Load ChEA3 results
tf_results <- read.table(
  "/Users/sowmi/Desktop/IPF_Research/ChEA3_TF_results.tsv",
  header = TRUE,
  sep = "\t"
)

# Check columns
head(tf_results)
colnames(tf_results)


# ---- CLEAN TF VISUALIZATION ----
library(ggplot2)

# Get top 15 TFs
tf_top15 <- head(tf_results, 15)

# Count overlapping genes for each TF
tf_top15$n_genes <- sapply(
  strsplit(tf_top15$Overlapping_Genes, ","), length)

# Sort for plotting
tf_top15 <- tf_top15[order(tf_top15$n_genes, decreasing = TRUE), ]
tf_top15$TF <- factor(tf_top15$TF, levels = rev(tf_top15$TF))

# Plot
ggplot(tf_top15, aes(x = TF, y = n_genes, fill = Score)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "firebrick", high = "steelblue",
                      name = "Mean Rank\n(lower = better)") +
  coord_flip() +
  labs(
    title    = "Top 15 Transcription Factors",
    subtitle = "Regulating 122 FUS-IPF Overlap Genes (ChEA3)",
    x = "Transcription Factor",
    y = "Number of Overlapping Genes"
  ) +
  theme_classic() +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    axis.text.y   = element_text(size = 11, face = "bold")
  )

ggsave("/Users/sowmi/Desktop/IPF_Research/TF_analysis.pdf",
       width = 9, height = 7)

# Save top TF results
write.csv(tf_top15,
          "/Users/sowmi/Desktop/IPF_Research/TF_top15_results.csv",
          row.names = FALSE)

cat("TF figure saved!\n")
cat("\nTop 5 TFs regulating your overlap genes:\n")
print(tf_top15[1:5, c("TF", "n_genes", "Score")])

# Install IOBR package for immune deconvolution
if (!require("IOBR", quietly = TRUE)) {
  if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
  devtools::install_github("IOBR/IOBR")
}

library(IOBR)

# Get normalized expression matrix from GSE99621
# We need TPM-like values — use VST normalized data
expr_matrix <- assay(vsd99)

# Check dimensions
dim(expr_matrix)
head(expr_matrix[1:5, 1:5])

library(IOBR)

# IOBR needs genes as rows, samples as columns
# Our matrix is already in this format — perfect

# Run CIBERSORT immune deconvolution
# This estimates 22 immune cell types in each sample
immune_results <- deconvo_tme(
  eset        = expr_matrix,
  method      = "cibersort",
  arrays      = FALSE,
  perm        = 100
)

# Check results
dim(immune_results)
head(immune_results[, 1:6])

library(ggplot2)
library(tidyr)
library(dplyr)

# Add condition labels to immune results
immune_df <- as.data.frame(immune_results)
immune_df$condition <- metadata99$condition

# Get immune cell columns only (remove ID and non-cell columns)
immune_cols <- grep("CIBERSORT", colnames(immune_df), value = TRUE)

# Remove p-value and correlation columns if present
immune_cols <- immune_cols[!grepl("p_value|correlation|RMSE", immune_cols)]

# Reshape to long format for plotting
immune_long <- immune_df %>%
  select(condition, all_of(immune_cols)) %>%
  pivot_longer(cols = all_of(immune_cols),
               names_to = "cell_type",
               values_to = "fraction") %>%
  mutate(cell_type = gsub("_CIBERSORT", "", cell_type),
         cell_type = gsub("_", " ", cell_type))

# Calculate mean per group
immune_summary <- immune_long %>%
  group_by(condition, cell_type) %>%
  summarise(mean_fraction = mean(fraction), .groups = "drop")

# Plot stacked bar chart
ggplot(immune_summary, 
       aes(x = condition, y = mean_fraction, fill = cell_type)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(
    title    = "Immune Cell Infiltration: IPF vs Control",
    subtitle = "CIBERSORT Deconvolution (GSE99621)",
    x        = "Condition",
    y        = "Mean Cell Fraction",
    fill     = "Cell Type"
  ) +
  theme_classic() +
  theme(
    plot.title    = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 11),
    legend.text   = element_text(size = 8),
    axis.text     = element_text(size = 12)
  )

ggsave("/Users/sowmi/Desktop/IPF_Research/immune_infiltration.pdf",
       width = 10, height = 8)

# Also find which cell types are significantly different
cat("\nMean immune fractions by condition:\n")
immune_summary %>%
  pivot_wider(names_from = condition, 
              values_from = mean_fraction) %>%
  mutate(difference = IPF - Control) %>%
  arrange(desc(abs(difference))) %>%
  print(n = 22)


#------------------------------------------------------------------------------

# Remove non-cell rows (Correlation, P-value, RMSE)
immune_clean <- immune_long %>%
  filter(!cell_type %in% c("Correlation", "P value", "RMSE"))

# Calculate statistics for each cell type
library(dplyr)

# Wilcoxon test for each cell type
cell_types_list <- unique(immune_clean$cell_type)

wilcox_results <- lapply(cell_types_list, function(ct) {
  ipf_vals  <- immune_clean$fraction[
    immune_clean$cell_type == ct & immune_clean$condition == "IPF"]
  ctrl_vals <- immune_clean$fraction[
    immune_clean$cell_type == ct & immune_clean$condition == "Control"]
  
  test <- wilcox.test(ipf_vals, ctrl_vals)
  data.frame(cell_type = ct, p_value = test$p.value)
})

wilcox_df <- do.call(rbind, wilcox_results)
wilcox_df$significance <- ifelse(wilcox_df$p_value < 0.05, "*", "")
wilcox_df <- wilcox_df[order(wilcox_df$p_value), ]

cat("Significantly different cell types (p < 0.05):\n")
print(wilcox_df[wilcox_df$p_value < 0.05, ])

# Boxplot of top significantly different cell types
sig_cells <- wilcox_df$cell_type[wilcox_df$p_value < 0.05]

if(length(sig_cells) > 0) {
  immune_sig <- immune_clean %>%
    filter(cell_type %in% sig_cells)
  
  ggplot(immune_sig, 
         aes(x = condition, y = fraction, fill = condition)) +
    geom_boxplot(outlier.shape = 21, alpha = 0.8) +
    geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
    facet_wrap(~ cell_type, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c("Control" = "steelblue",
                                 "IPF"     = "firebrick")) +
    labs(
      title    = "Significantly Different Immune Cell Types",
      subtitle = "IPF vs Control — CIBERSORT (Wilcoxon p < 0.05)",
      x = "", y = "Cell Fraction"
    ) +
    theme_classic() +
    theme(
      strip.text    = element_text(size = 9, face = "bold"),
      plot.title    = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.position = "none"
    )
  
  ggsave("/Users/sowmi/Desktop/IPF_Research/immune_boxplots.pdf",
         width = 12, height = 8)
  cat("Immune boxplot saved!\n")
} else {
  cat("No significant differences at p<0.05 — try p<0.1\n")
}


# Remove non-cell type rows and replot
immune_clean2 <- immune_long %>%
  filter(!cell_type %in% c("Correlation", "P value", 
                           "RMSE", "P-value"))

# Redo wilcoxon on clean data
cell_types_list2 <- unique(immune_clean2$cell_type)

wilcox_results2 <- lapply(cell_types_list2, function(ct) {
  ipf_vals  <- immune_clean2$fraction[
    immune_clean2$cell_type == ct & 
      immune_clean2$condition == "IPF"]
  ctrl_vals <- immune_clean2$fraction[
    immune_clean2$cell_type == ct & 
      immune_clean2$condition == "Control"]
  test <- wilcox.test(ipf_vals, ctrl_vals)
  data.frame(cell_type = ct, p_value = test$p.value)
})

wilcox_df2 <- do.call(rbind, wilcox_results2)
sig_cells2 <- wilcox_df2$cell_type[wilcox_df2$p_value < 0.05]

cat("Significant cell types:", paste(sig_cells2, collapse=", "), "\n")

# Clean boxplot
immune_sig2 <- immune_clean2 %>%
  filter(cell_type %in% sig_cells2)

ggplot(immune_sig2,
       aes(x = condition, y = fraction, fill = condition)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.8) +
  geom_jitter(width = 0.2, size = 1.5, alpha = 0.6) +
  facet_wrap(~ cell_type, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = c("Control" = "steelblue",
                               "IPF"     = "firebrick")) +
  labs(
    title    = "Significantly Different Immune Cell Types",
    subtitle = "IPF vs Control — CIBERSORT (Wilcoxon p < 0.05)",
    x = "", y = "Cell Fraction"
  ) +
  theme_classic() +
  theme(
    strip.text      = element_text(size = 10, face = "bold"),
    plot.title      = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle   = element_text(hjust = 0.5, size = 11),
    legend.position = "none"
  )

ggsave("/Users/sowmi/Desktop/IPF_Research/immune_boxplots_clean.pdf",
       width = 10, height = 6)

# Save immune results
write.csv(as.data.frame(immune_results),
          "/Users/sowmi/Desktop/IPF_Research/CIBERSORT_results.csv",
          row.names = FALSE)

cat("Clean immune figures saved!\n")



# ---- PREPARE DATA FOR ML ----
library(DESeq2)

# Get VST normalized expression for GSE99621
# Keep only the 122 overlap genes
vsd99_matrix <- assay(vsd99)

# Find which overlap genes are in the expression matrix
overlap_in_matrix <- overlap_genes[
  overlap_genes %in% rownames(vsd99_matrix)]

cat("Overlap genes found in matrix:", 
    length(overlap_in_matrix), "\n")

# Extract those genes
overlap_expr <- vsd99_matrix[overlap_in_matrix, ]

# Transpose — ML needs samples as rows, genes as columns
ml_data <- as.data.frame(t(overlap_expr))

# Add condition as binary outcome (1=IPF, 0=Control)
ml_data$condition <- ifelse(
  metadata99$condition == "IPF", 1, 0)

# Check
cat("Samples:", nrow(ml_data), "\n")
cat("Features:", ncol(ml_data) - 1, "\n")
cat("IPF:", sum(ml_data$condition == 1), "\n")
cat("Control:", sum(ml_data$condition == 0), "\n")

# ---- STEP 1: LASSO REGRESSION ----
install.packages("glmnet")
library(glmnet)

set.seed(42)

# Prepare X (features) and y (outcome)
X <- as.matrix(ml_data[, -ncol(ml_data)])
y <- ml_data$condition

# Run LASSO with cross-validation
lasso_cv <- cv.glmnet(
  x          = X,
  y          = y,
  family     = "binomial",
  alpha      = 1,        # 1 = LASSO
  nfolds     = 5,        # 5-fold cross validation
  standardize = TRUE
)

# Plot cross-validation curve
plot(lasso_cv)
title("LASSO Cross-Validation", line = 3)

# Get best lambda
best_lambda <- lasso_cv$lambda.min
cat("Best lambda:", best_lambda, "\n")

# Get selected genes at best lambda
lasso_coef <- coef(lasso_cv, s = "lambda.min")
selected_genes <- rownames(lasso_coef)[
  which(abs(lasso_coef) > 0)]
selected_genes <- selected_genes[
  selected_genes != "(Intercept)"]

cat("Genes selected by LASSO:", length(selected_genes), "\n")
cat("Selected genes:\n")
print(selected_genes)


# ---- STEP 2: RANDOM FOREST + SHAP ----
install.packages(c("randomForest", "shapviz", "kernelshap"))
library(randomForest)

set.seed(42)

# Train Random Forest on the 7 LASSO-selected genes
rf_data <- ml_data[, c(selected_genes, "condition")]
rf_data$condition <- as.factor(rf_data$condition)

# Train model
rf_model <- randomForest(
  condition ~ .,
  data       = rf_data,
  ntree      = 500,
  importance = TRUE
)

print(rf_model)

# Variable importance plot
varImpPlot(rf_model,
           main = "Random Forest Feature Importance\n7 LASSO-Selected Genes")

# Get importance values
importance_df <- as.data.frame(importance(rf_model))
importance_df$gene <- rownames(importance_df)
importance_df <- importance_df[order(
  importance_df$MeanDecreaseGini, decreasing=TRUE), ]

cat("\nGene Importance Ranking:\n")
print(importance_df[, c("gene", "MeanDecreaseGini")])



# ---- STEP 3: SVM WITH 5-FOLD CROSS VALIDATION ----
install.packages(c("e1071", "pROC", "caret"))
library(e1071)
library(pROC)
library(caret)

set.seed(42)

# Prepare data with selected genes only
svm_data <- ml_data[, c(selected_genes, "condition")]
svm_data$condition <- as.factor(svm_data$condition)

# 5-fold cross validation setup
folds <- createFolds(svm_data$condition, k = 5, list = TRUE)

# Store predictions
all_probs  <- c()
all_labels <- c()

for(i in 1:5) {
  # Split train/test
  test_idx  <- folds[[i]]
  train_idx <- unlist(folds[-i])
  
  train_data <- svm_data[train_idx, ]
  test_data  <- svm_data[test_idx, ]
  
  # Train SVM
  svm_model <- svm(
    condition ~ .,
    data        = train_data,
    kernel      = "radial",
    probability = TRUE,
    cost        = 1
  )
  
  # Predict
  pred <- predict(svm_model, test_data, probability = TRUE)
  probs <- attr(pred, "probabilities")[, "1"]
  
  all_probs  <- c(all_probs, probs)
  all_labels <- c(all_labels, 
                  as.numeric(as.character(test_data$condition)))
}

# Calculate ROC and AUC
roc_obj <- roc(all_labels, all_probs)
auc_val <- auc(roc_obj)
cat("AUC:", round(auc_val, 4), "\n")

# Plot ROC curve
pdf("/Users/sowmi/Desktop/IPF_Research/ROC_curve.pdf",
    width = 7, height = 7)
plot(roc_obj,
     col  = "firebrick",
     lwd  = 2,
     main = paste0("ROC Curve — 7-Gene IPF Signature\n",
                   "AUC = ", round(auc_val, 4)),
     print.auc = TRUE)
abline(a = 0, b = 1, lty = 2, col = "grey50")
dev.off()

cat("ROC curve saved!\n")

set.seed(42)

library(glmnet)
library(randomForest)
library(pROC)
library(caret)

# Prepare data
X <- as.matrix(ml_data[, -ncol(ml_data)])
y <- as.factor(ml_data$condition)

# 5-fold cross validation — PROPER version
folds <- createFolds(y, k = 5, list = TRUE)

all_probs  <- c()
all_labels <- c()
all_selected_genes <- list()

for(i in 1:5) {
  
  cat("Fold", i, "...\n")
  
  # Split FIRST
  test_idx  <- folds[[i]]
  train_idx <- unlist(folds[-i])
  
  X_train <- X[train_idx, ]
  X_test  <- X[test_idx, ]
  y_train <- y[train_idx]
  y_test  <- y[test_idx]
  
  # LASSO on TRAINING data only
  lasso_fold <- cv.glmnet(
    x       = X_train,
    y       = as.numeric(as.character(y_train)),
    family  = "binomial",
    alpha   = 1,
    nfolds  = 3
  )
  
  # Get genes selected in this fold
  coef_fold <- coef(lasso_fold, s = "lambda.min")
  genes_fold <- rownames(coef_fold)[
    which(abs(coef_fold) > 0)]
  genes_fold <- genes_fold[genes_fold != "(Intercept)"]
  all_selected_genes[[i]] <- genes_fold
  
  if(length(genes_fold) == 0) {
    cat("No genes selected in fold", i, "- skipping\n")
    next
  }
  
  # Random Forest on TRAINING data with selected genes
  train_df <- as.data.frame(X_train[, genes_fold, drop=FALSE])
  train_df$condition <- y_train
  
  test_df  <- as.data.frame(X_test[, genes_fold, drop=FALSE])
  
  rf_fold <- randomForest(
    condition ~ .,
    data  = train_df,
    ntree = 500,
    probability = TRUE
  )
  
  # Predict on TEST data (never seen by model)
  pred_fold <- predict(rf_fold, test_df, type = "prob")
  probs_fold <- pred_fold[, "1"]
  
  all_probs  <- c(all_probs, probs_fold)
  all_labels <- c(all_labels, 
                  as.numeric(as.character(y_test)))
}

# Real AUC
roc_obj <- roc(all_labels, all_probs)
auc_val <- auc(roc_obj)
cat("\nReal AUC (proper cross-validation):", 
    round(auc_val, 4), "\n")

# Which genes were consistently selected?
cat("\nGenes selected across folds:\n")
gene_freq <- table(unlist(all_selected_genes))
gene_freq <- sort(gene_freq, decreasing = TRUE)
print(gene_freq)
cat("\nGenes selected in ALL 5 folds:\n")
print(names(gene_freq[gene_freq == 5]))

# Fix gene names with hyphens — replace with underscores
colnames(ml_data) <- make.names(colnames(ml_data))

# Rebuild X and y with fixed names
X <- as.matrix(ml_data[, -ncol(ml_data)])
y <- as.factor(ml_data$condition)

set.seed(42)
folds <- createFolds(y, k = 5, list = TRUE)

all_probs  <- c()
all_labels <- c()
all_selected_genes <- list()

for(i in 1:5) {
  cat("Fold", i, "...\n")
  
  test_idx  <- folds[[i]]
  train_idx <- unlist(folds[-i])
  
  X_train <- X[train_idx, ]
  X_test  <- X[test_idx, ]
  y_train <- y[train_idx]
  y_test  <- y[test_idx]
  
  # LASSO on training only
  lasso_fold <- cv.glmnet(
    x      = X_train,
    y      = as.numeric(as.character(y_train)),
    family = "binomial",
    alpha  = 1,
    nfolds = 3
  )
  
  coef_fold  <- coef(lasso_fold, s = "lambda.min")
  genes_fold <- rownames(coef_fold)[which(abs(coef_fold) > 0)]
  genes_fold <- genes_fold[genes_fold != "(Intercept)"]
  all_selected_genes[[i]] <- genes_fold
  
  if(length(genes_fold) == 0) {
    cat("No genes in fold", i, "- skipping\n")
    next
  }
  
  # Random Forest on training only
  train_df <- as.data.frame(X_train[, genes_fold, drop=FALSE])
  train_df$condition <- y_train
  test_df  <- as.data.frame(X_test[, genes_fold, drop=FALSE])
  
  rf_fold <- randomForest(
    condition ~ .,
    data  = train_df,
    ntree = 500
  )
  
  pred_fold  <- predict(rf_fold, test_df, type = "prob")
  probs_fold <- pred_fold[, "1"]
  
  all_probs  <- c(all_probs, probs_fold)
  all_labels <- c(all_labels, as.numeric(as.character(y_test)))
}

# Real AUC
roc_obj <- roc(all_labels, all_probs)
auc_val <- auc(roc_obj)
cat("\nHonest AUC (proper cross-validation):", round(auc_val, 4), "\n")

# Genes consistently selected
cat("\nGenes selected per fold:\n")
gene_freq <- sort(table(unlist(all_selected_genes)), decreasing=TRUE)
print(gene_freq)
cat("\nGenes selected in 3+ folds (most stable):\n