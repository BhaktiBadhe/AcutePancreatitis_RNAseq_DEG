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

# ----- Visualization
theme_bold <- function(base_size = 16) {
  theme_classic(base_size = base_size) +
    theme(
      axis.title.x = element_text(face = "bold", size = base_size, color = "black"),
      axis.title.y = element_text(face = "bold", size = base_size, color = "black"),
      axis.text.x  = element_text(face = "bold", size = base_size - 2, color = "black"),
      axis.text.y  = element_text(face = "bold", size = base_size - 2, color = "black"),
      
      plot.title   = element_text(face = "bold", size = base_size + 2, hjust = 0.5),
      
      legend.title = element_text(face = "bold", size = base_size - 2),
      legend.text  = element_text(size = base_size - 3, face = "bold"),
      legend.position = "right")
}

# PCA plot
# Variance stabilizing transformation
vsd <- vst(dds, blind = FALSE)

# PCA
ggplot(pca_df, aes(PC1, PC2, color = group2)) +
  geom_point(size = 4) +
  xlab(paste0("PC1 (", percentVar[1], "% variance)")) +
  ylab(paste0("PC2 (", percentVar[2], "% variance)")) +
  ggtitle("Acute Pancreatitis vs Healthy Controls") +
  guides(color = guide_legend(title = "Group")) +
  theme_bold()

# Volcano plot
volcano_df <- res_deg_df %>%
  mutate(sig = case_when(
      padj < 0.05 & log2FoldChange > 1  ~ "Upregulated",
      padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
      TRUE ~ "Not significant"))

ggplot(volcano_df, aes(log2FoldChange, -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c(
    "Upregulated" = "red",
    "Downregulated" = "blue",
    "Not significant" = "grey")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic() +
  labs(title = "AP vs HC",
       x = "log2 Fold Change",
       y = "-log10 adjusted P-value") + theme_bold()

# Heatmap for top DEGs
# Select top 10 DEGs
top_genes <- sig_deg %>%
  slice_head(n = 10) %>%
  pull(gene)

# Extract normalized expression
mat <- assay(vsd)[top_genes, ]

# Scale by gene
mat_scaled <- t(scale(t(mat)))

# Annotation
annotation_col <- data.frame(Group = colData(dds)$group2)
rownames(annotation_col) <- colnames(mat_scaled)

# Heatmap
pheatmap( mat_scaled,
  annotation_col = annotation_col,
  show_colnames = FALSE,
  fontsize_row = 8,
  main = "Top 10 DEGs") + theme_bold()









