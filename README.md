# 🧬 Integrative Analysis of FUS-Regulated Gene Networks in Idiopathic Pulmonary Fibrosis (IPF)

## 📌 Overview

This project investigates the role of the RNA-binding protein **FUS** in regulating gene expression associated with **Idiopathic Pulmonary Fibrosis (IPF)**.

By integrating **perturbation-based RNA-seq data (FUS knockdown)** with **disease-level transcriptomic data**, this study identifies key gene networks, pathways, and regulatory mechanisms involved in fibrotic progression.

---

## 🎯 Objectives

* Identify genes regulated by FUS using RNA-seq data
* Detect differentially expressed genes (DEGs) in IPF
* Perform cross-dataset validation to identify **common genes**
* Analyze functional pathways (GO, KEGG)
* Construct protein-protein interaction (PPI) networks
* Identify regulatory elements (Transcription Factors and miRNAs)
* Explore potential spatial organization using Hi-C data (ongoing)

---

## 🧪 Datasets Used

* **GSE310475** – FUS knockdown (ASO vs Control)
* **E-GEOD-52463** – IPF vs healthy lung tissue

---

## ⚙️ Methodology

### 1. Data Processing

* Raw count data loaded and cleaned in R
* Sample grouping: Control vs FUS knockdown / IPF vs healthy

### 2. Differential Expression Analysis

* Performed using **DESeq2**
* Thresholds:

  * Adjusted p-value < 0.05
  * |log2FoldChange| > 1

### 3. Cross-Dataset Integration

* Intersection of DEGs from both datasets
* Identification of **122 overlapping genes**

### 4. Functional Enrichment

* Gene Ontology (GO) analysis
* KEGG pathway analysis

### 5. Network Analysis

* Protein-Protein Interaction (PPI) networks
* Hub gene identification

### 6. Regulatory Analysis

* Transcription Factor prediction (e.g., TCF21, TP63)
* miRNA network analysis (e.g., miR-34a, let-7b)

---

## 📊 Key Results

* **633 DEGs** identified in FUS knockdown dataset
* **1486 DEGs** identified in IPF dataset
* **122 overlapping genes** between FUS and IPF

### 🔬 Biological Insights

* Enriched pathways include:

  * Extracellular Matrix (ECM) organization
  * TGF-β signaling
  * Immune response pathways

### 🧬 Regulatory Insights

* Key transcription factors: **TCF21, TP63**
* Key miRNAs: **miR-34a, let-7b**

---

## 🧠 Key Conclusion

The results suggest that **FUS plays a regulatory role in IPF progression** by modulating gene networks involved in fibrosis, immune response, and extracellular matrix remodeling.

---

## 🔭 Future Work

* Integration of **Hi-C data** to study 3D genome organization
* Validation using additional datasets
* Experimental validation of key genes and pathways

---

## 🛠️ Tools & Technologies

* R
* DESeq2
* clusterProfiler
* STRING database
* ggplot2 / visualization tools

---

## 📁 Project Structure

```
IPF_Research/
│── data/               # Raw and processed datasets
│── scripts/            # R scripts for analysis
│── results/            # DEG lists, enrichment results
│── figures/            # Plots and visualizations
│── README.md           # Project documentation
```

---

## 📌 Author

**Sowmya Arunachalam**
Integrated M.Tech (CSE – Business Analytics)
Interested in AI, Bioinformatics, and Computational Biology

**Cavin Chandran CA**
Integrated M.Tech (CSE – Business Analytics)

---

## ⭐ Acknowledgment

This project was developed as part of an independent exploration into transcriptomics and disease modeling, with a focus on applying computational techniques to real-world biological problems.

---
