# Differential Gene Expression Analysis of Acute Pancreatitis (AP) vs Healthy Controls

This repository contains an RNA-seq differential gene expression analysis comparing **Acute Pancreatitis (AP)** samples (moderate + severe) with **Healthy Controls (HC)** using the **DESeq2** framework in R.

The analysis is based on publicly available count data from **GEO accession GSE194331**.

---

## 📌 Study Design

- **Organism:** Homo sapiens  
- **Conditions analyzed:**
  - Healthy Controls (HC)
  - Acute Pancreatitis (AP: Moderate + Severe)
- **Mild AP samples were excluded**
- **All samples originate from the same study**, minimizing batch effects

---

## ⚙️ Methods Overview

### 1. Data Input
- Raw gene-level count matrix
- Sample metadata containing clinical condition labels

### 2. Pre-processing
- Standardization of condition labels
- Removal of mild AP samples
- Filtering of lowly expressed genes (row sum > 10)

### 3. Differential Expression Analysis
- DESeq2 size-factor normalization (median-of-ratios method)
- Negative binomial generalized linear model
- Contrast tested: **AP vs HC**
- Significance thresholds:
  - Adjusted p-value < 0.05
  - |log2 fold change| > 1

### 4. Visualization
- Variance Stabilizing Transformation (VST)
- Principal Component Analysis (PCA)
- Volcano plot
- Heatmap of top 10 differentially expressed genes

---

## Output Files

| File | Description |
|-----|-------------|
| `AP_vs_HC_DEGs.csv` | Significantly differentially expressed genes |
| `PCA_plot.png` | PCA plot using VST-transformed data |
| `Volcano_plot.png` | Volcano plot of AP vs HC |
| `Heatmap.png` | Heatmap of top 10 DEGs |

---

## 🧬 R Packages Used

- DESeq2
- tidyverse
- pheatmap
- stringr
- ggplot2

---

## Run the Analysis

```r
# Install required packages
install.packages(c("tidyverse", "pheatmap", "stringr"))
BiocManager::install("DESeq2")

# Run the analysis
source("AP_deg_pipeline.R")

