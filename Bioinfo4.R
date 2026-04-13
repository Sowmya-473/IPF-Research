# Check what's still in memory
ls()


# ---- GSEA on GSE310475 (FUS silencing dataset) ----
install.packages("msigdbr")
BiocManager::install("fgsea")

library(fgsea)
library(msigdbr)
library(dplyr)
library(ggplot2)

# Step 1: Create ranked gene list (all genes, not just significant ones)
# Rank by stat (more reliable than log2FC alone)
ranked_gse <- res_df$stat
names(ranked_gse) <- res_df$gene_symbol

# Remove NAs and duplicates
ranked_gse <- ranked_gse[!is.na(names(ranked_gse))]
ranked_gse <- ranked_gse[!duplicated(names(ranked_gse))]

# Sort descending
ranked_gse <- sort(ranked_gse, decreasing = TRUE)

cat("Ranked gene list length:", length(ranked_gse), "\n")
cat("Top 5 genes:", paste(names(head(ranked_gse, 5)), collapse=", "), "\n")
cat("Bottom 5 genes:", paste(names(tail(ranked_gse, 5)), collapse=", "), "\n")


# Step 2: Get gene sets from MSigDB
# We'll use three collections:
# H  = Hallmark (50 curated pathways, most interpretable)
# C2 = KEGG/Reactome (canonical pathways)
# C5 = GO biological process

# Fixed version for msigdbr 10.0.0+
h_sets <- msigdbr(species = "Homo sapiens", 
                  collection = "H") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Correct names for msigdbr 2026
c2_sets <- msigdbr(species = "Homo sapiens",
                   collection = "C2",
                   subcollection = "CP:KEGG_LEGACY") %>%
  split(x = .$gene_symbol, f = .$gs_name)

c5_sets <- msigdbr(species = "Homo sapiens",
                   collection = "C5",
                   subcollection = "GO:BP") %>%
  split(x = .$gene_symbol, f = .$gs_name)

# Also grab Reactome — more comprehensive than KEGG
reactome_sets <- msigdbr(species = "Homo sapiens",
                         collection = "C2",
                         subcollection = "CP:REACTOME") %>%
  split(x = .$gene_symbol, f = .$gs_name)

cat("Hallmark gene sets:", length(h_sets), "\n")
cat("KEGG gene sets:", length(c2_sets), "\n")
cat("Reactome gene sets:", length(reactome_sets), "\n")
cat("GO:BP gene sets:", length(c5_sets), "\n")

# Now run GSEA on all four
set.seed(42)

gsea_hallmark <- fgsea(pathways = h_sets,
                       stats    = ranked_gse,
                       minSize  = 15,
                       maxSize  = 500,
                       nperm    = 1000)

gsea_kegg <- fgsea(pathways = c2_sets,
                   stats    = ranked_gse,
                   minSize  = 15,
                   maxSize  = 500,
                   nperm    = 1000)

gsea_reactome <- fgsea(pathways = reactome_sets,
                       stats    = ranked_gse,
                       minSize  = 15,
                       maxSize  = 500,
                       nperm    = 1000)

gsea_go <- fgsea(pathways = c5_sets,
                 stats    = ranked_gse,
                 minSize  = 15,
                 maxSize  = 500,
                 nperm    = 1000)

# Filter significant results
sig_hallmark  <- gsea_hallmark  %>% filter(padj < 0.05) %>% arrange(NES)
sig_kegg      <- gsea_kegg      %>% filter(padj < 0.05) %>% arrange(NES)
sig_reactome  <- gsea_reactome  %>% filter(padj < 0.05) %>% arrange(NES)
sig_go        <- gsea_go        %>% filter(padj < 0.05) %>% arrange(NES)

cat("\nSignificant Hallmark pathways:", nrow(sig_hallmark), "\n")
cat("Significant KEGG pathways:", nrow(sig_kegg), "\n")
cat("Significant Reactome pathways:", nrow(sig_reactome), "\n")
cat("Significant GO:BP terms:", nrow(sig_go), "\n")

# Show Hallmark results — most interpretable
cat("\n--- HALLMARK pathways (NES < 0 = suppressed by FUS silencing) ---\n")
print(sig_hallmark[, c("pathway", "NES", "padj")])