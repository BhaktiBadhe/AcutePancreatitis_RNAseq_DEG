# set the path to the folder
# ------------------------------------------------------------------------------------

library(DESeq2)
library(tidyverse)
library(pheatmap)
library(stringr)

# Load counts
counts <- read.table("GSE194331_HC_PAN_PANSEP_counts.txt",header = TRUE,
  row.names = 1, sep = "\t", check.names = FALSE)

# Load metadata
metadata <- read.csv("metadata.csv", header = TRUE, stringsAsFactors = FALSE)

# Standardize condition labels
metadata$condition <- case_when(
  str_trim(metadata$Condition) == "Healthy control" ~ "HC",
  str_trim(metadata$Condition) == "Moderately-severe AP" ~ "Moderate",
  str_trim(metadata$Condition) == "Severe AP" ~ "Severe",
  str_trim(metadata$Condition) == "Mild AP" ~ "Mild")

# Drop mild AP
metadata <- metadata %>%
  filter(condition %in% c("HC", "Moderate", "Severe"))

# Match samples
common_samples <- intersect(colnames(counts), metadata$Sample)

counts <- counts[, common_samples]
metadata <- metadata[match(common_samples, metadata$Sample), ]

# verify order
stopifnot(all(colnames(counts) == metadata$Sample))

# Finalize metadata
rownames(metadata) <- metadata$Sample

metadata$condition <- factor(metadata$condition,
  levels = c("HC", "Moderate", "Severe"))
table(metadata$condition)

# Create DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = counts,
  colData = metadata,
  design = ~ condition)   # dim: 58735    62

# Filter low count genes
dds <- dds[rowSums(counts(dds)) > 10, ]  # dim: 26915 62

# Collapse Moderate + Severe into AP
colData(dds)$group2 <- ifelse(colData(dds)$condition == "HC", "HC", "AP")
colData(dds)$group2 <- factor(colData(dds)$group2, levels = c("HC", "AP"))
table(colData(dds)$group2)           

# Update design and run DESeq2
design(dds) <- ~ group2
dds <- DESeq(dds)

# Differential expression: AP vs HC
res_deg <- results(dds, contrast = c("group2", "AP", "HC"))
summary(res_deg)

# Clean DEG table
res_deg_df <- as.data.frame(res_deg) %>%
  rownames_to_column("gene") %>%
  arrange(padj)

sig_deg <- res_deg_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1)

nrow(sig_deg)

# Save results
write.csv(sig_deg,"AP_vs_HC_DEGs.csv",row.names = FALSE)


