# ===============================
# 1. Set working directory (Mac)
# ===============================
setwd("/Users/brindha/Downloads/Mini Project/M1")

# ===============================
# 2. Load count data
# ===============================
counts <- read.table("GSE164416_DP_htseq_counts.txt",
                     header = TRUE,
                     row.names = 1,
                     sep = "\t",
                     check.names = FALSE)

# Quick checks
cat("Counts dimensions: ", dim(counts), "\n")

# ===============================
# 3. Parse metadata from GEO file
# ===============================
file_name <- "GSE164416_series_matrix.txt"
lines <- readLines(file_name)

# Extract GSM IDs
gsm_line <- grep("^!Sample_geo_accession", lines, value = TRUE)
gsm_ids <- unlist(strsplit(gsm_line, "\t"))[-1]

# Extract sample titles
title_line <- grep("^!Sample_title", lines, value = TRUE)
titles <- unlist(strsplit(title_line, "\t"))[-1]

# Extract DP IDs and condition labels
library(stringr)
dp_ids <- str_extract(titles, "DP\\d+")
conditions <- str_extract(titles, "T3cD|T2D|IGT|ND|IFG")

# Build metadata
metadata <- data.frame(
  sample = dp_ids,
  GSM = gsm_ids,
  condition = conditions,
  stringsAsFactors = FALSE
)

# Reorder to match count matrix
metadata <- metadata[match(colnames(counts), metadata$sample), ]
rownames(metadata) <- metadata$sample

# ===============================
# 4. Save metadata to Excel
# ===============================
if (!requireNamespace("openxlsx", quietly = TRUE)) {
  install.packages("openxlsx")
}
library(openxlsx)

write.xlsx(metadata, "Cleaned_Metadata.xlsx", rowNames = TRUE)

cat("✅ Cleaned metadata saved as: Cleaned_Metadata.xlsx in working directory\n")

# ===============================
# 5. Check unique conditions
# ===============================
print(unique(metadata$condition))

# ===============================
# 6. Load libraries
# ===============================
if (!requireNamespace("DESeq2", quietly = TRUE)) {
  install.packages("BiocManager")
  BiocManager::install("DESeq2")
}
if (!requireNamespace("pheatmap", quietly = TRUE)) {
  install.packages("pheatmap")
}
if (!requireNamespace("RColorBrewer", quietly = TRUE)) {
  install.packages("RColorBrewer")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  install.packages("ggplot2")
}

library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

# ===============================
# 7. Create DESeq2 object
# ===============================
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ condition)

# Filter out low-count genes
dds <- dds[rowSums(counts(dds)) > 10, ]

# Variance stabilizing transformation
vsd <- vst(dds, blind = TRUE)

# ===============================
# 8. Heatmap of sample distances
# ===============================
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$condition
colnames(sampleDistMatrix) <- vsd$condition

pheatmap(sampleDistMatrix,
         clustering_distance_rows = sampleDists,
         clustering_distance_cols = sampleDists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-Sample Distance")

# ===============================
# 9. Heatmap of top variable genes
# ===============================
topVarGenes <- head(order(rowVars(assay(vsd)), decreasing = TRUE), 50)
pheatmap(assay(vsd)[topVarGenes, ],
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         annotation_col = metadata,
         col = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
         main = "Top 50 Variable Genes")

# ===============================
# 10. PCA Plot
# ===============================
pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of RNA-seq Samples") +
  theme_minimal()

# ===============================
# 11. Clustering Dendrogram
# ===============================
hc <- hclust(sampleDists)
plot(hc, main = "Hierarchical Clustering of Samples", xlab = "", sub = "")



# ===============================
# 12. Box plot
# ===============================
boxplot(assay(vsd),
        col = RColorBrewer::brewer.pal(8, "Set3"),
        main = "Boxplot of Normalized Counts",
        las = 2, outline = FALSE,
        ylab = "VST normalized counts")

# ===============================
# 12. Density plot of expression distributions
# ===============================
# log-transform counts
log_counts <- log2(counts + 1)

# initialize plot with first sample
d <- density(log_counts[,1])
plot(d$x, d$y, type="l", col=1, lwd=2,
     main="Density of log2(counts+1)",
     xlab="log2(Expression+1)", ylab="Density")

# add remaining samples
for(i in 2:ncol(log_counts)) {
  d <- density(log_counts[,i])
  lines(d$x, d$y, col=i, lwd=2)
}

# add legend
legend("topright", legend=colnames(log_counts),
       col=1:ncol(log_counts), lty=1, cex=0.7)

# ===============================
# 12. biomarker discovery
# ===============================
# ===============================
# BIOMARKER SCREEN: T2D vs ND
# Base R + optional pROC (CRAN)
# ===============================

# 0) Safety checks and alignments
stopifnot(all(colnames(counts) %in% rownames(metadata)))
metadata <- metadata[colnames(counts), , drop = FALSE]

# 1) Choose comparison groups
groupA <- "T2D"  # case
groupB <- "ND"   # control

keep <- metadata$condition %in% c(groupA, groupB)
if (!any(keep)) stop("No samples for the chosen groups in metadata$condition.")
counts_sub <- counts[, keep, drop = FALSE]
md_sub     <- metadata[keep, , drop = FALSE]
group      <- factor(ifelse(md_sub$condition == groupA, groupA, groupB), levels = c(groupB, groupA))

# 2) Transform counts (log2 CPM-like, simple and robust)
log_counts <- log2(counts_sub + 1)

# 3) Per-gene Wilcoxon test, log2FC, AUC (derived from Wilcoxon U)
nA <- sum(group == groupA)
nB <- sum(group == groupB)
if (nA < 2 || nB < 2) warning("One of the groups has <2 samples; stats may be unstable.")

wilcoxon_stats <- function(x, g) {
  xa <- x[g == groupA]
  xb <- x[g == groupB]
  # Wilcoxon rank-sum
  wt <- suppressWarnings(wilcox.test(xa, xb, exact = FALSE))
  # Wilcoxon statistic (W) is sum of ranks for first vector (xa)
  # Convert to U: U = W - nA*(nA+1)/2
  # AUC = U / (nA*nB)
  # Retrieve W from wt$statistic
  W  <- as.numeric(wt$statistic)
  U  <- W - nA * (nA + 1) / 2
  auc <- U / (nA * nB)
  # log2 FC (A - B)
  lfc <- mean(xa, na.rm = TRUE) - mean(xb, na.rm = TRUE)
  c(p = wt$p.value, log2FC = lfc, AUC = auc,
    
    # ===============================
    # 7) PCA (ggplot2) on top variable genes
    # ===============================
    
    # Select top variable genes (e.g., top 500)
    var_genes <- head(order(apply(log_counts, 1, var), decreasing = TRUE), 500)
    expr_top  <- log_counts[var_genes, , drop = FALSE]
    
    # Run PCA (samples as rows)
    pca <- prcomp(t(expr_top), scale. = TRUE)
    
    # Extract scores
    pca_df <- data.frame(pca$x, condition = group)
    
    # Calculate explained variance
    percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
    
    # Create color palette
    cols <- c("ND" = "#3B82F6", "T2D" = "#EF4444")
    
    # Plot
    library(ggplot2)
    png("PCA_T2D_vs_ND_colored.png", width = 1400, height = 1000, res = 130)
    ggplot(pca_df, aes(PC1, PC2, color = condition)) +
      geom_point(size = 5, alpha = 0.8) +
      xlab(paste0("PC1 (", percentVar[1], "% variance)")) +
      ylab(paste0("PC2 (", percentVar[2], "% variance)")) +
      ggtitle("PCA of Top 500 Variable Genes (T2D vs ND)") +
      scale_color_manual(values = cols) +
      theme_minimal(base_size = 18) +
      theme(
        legend.position = "top",
        plot.title = element_text(hjust = 0.5, face = "bold"),
        axis.title = element_text(face = "bold"),
        panel.grid = element_line(color = "gray85")
      ) +
      geom_text(aes(label = rownames(pca_df)), vjust = -1.2, size = 4, show.legend = FALSE)
    dev.off()
    cat("✅ Saved: PCA_T2D_vs_ND_colored.png (ggplot2 enhanced)\n")
    

    meanA = mean(xa, na.rm = TRUE), meanB = mean(xb, na.rm = TRUE))
}

res_mat <- t(apply(log_counts, 1, wilcoxon_stats, g = group))
res_df <- data.frame(
  gene   = rownames(log_counts),
  pvalue = res_mat[, "p"],
  log2FC = res_mat[, "log2FC"],
  AUC    = res_mat[, "AUC"],
  mean_case = res_mat[, "meanA"],
  mean_ctrl = res_mat[, "meanB"],
  stringsAsFactors = FALSE
)
res_df$padj <- p.adjust(res_df$pvalue, method = "BH")

# 4) Rank and save results
res_df <- res_df[order(res_df$padj, -abs(res_df$log2FC)), ]
write.csv(res_df, "Biomarker_screen_T2D_vs_ND.csv", row.names = FALSE)
cat("✅ Saved: Biomarker_screen_T2D_vs_ND.csv\n")

# 5) Volcano plot (base R)
sig <- with(res_df, padj < 0.05 & abs(log2FC) > 1)
col_vec <- ifelse(res_df$padj < 0.05 & res_df$log2FC > 1, "red",
                  ifelse(res_df$padj < 0.05 & res_df$log2FC < -1, "blue", "grey"))

png("Volcano_T2D_vs_ND_colored.png", width = 1200, height = 900, res = 130)
plot(res_df$log2FC, -log10(pmax(res_df$pvalue, .Machine$double.xmin)),
     xlab = "log2 Fold Change (T2D vs ND)",
     ylab = "-log10(p-value)",
     main = "Volcano: T2D vs ND",
     pch = 20, col = col_vec)
abline(v = c(-1, 1), lty = 2)
abline(h = -log10(0.05), lty = 2)
legend("topright", legend = c("Upregulated in T2D", "Downregulated in T2D", "Not sig."),
       col = c("red", "blue", "grey"), pch = 20, bty = "n")
dev.off()
cat("✅ Saved: Volcano_T2D_vs_ND_colored.png\n")

# 6) Heatmap of top genes (base heatmap)
####1 
head(res_df)
####2
topN <- 50
top_genes <- head(res_df$gene, topN)

if (length(top_genes) == 0) {
  stop("⚠️ No genes available for heatmap. Check your res_df content.")
}

mat_top <- log_counts[top_genes, , drop = FALSE]

# scale rows for better heatmap contrast
mat_top <- scale(t(scale(t(mat_top))))

ord_cols <- order(group)

####3
cols_group <- ifelse(group == groupA, "red", "blue")

heatmap(as.matrix(mat_top[, ord_cols]),
        Rowv = TRUE, Colv = NA, scale = "none",
        labCol = paste(colnames(mat_top)[ord_cols], as.character(group)[ord_cols], sep = "_"),
        col = colorRampPalette(c("navy", "white", "firebrick3"))(50),
        ColSideColors = cols_group[ord_cols],
        main = "Top 50 candidate biomarkers (scaled)")



# 7) PCA (base R) on top variable genes
cols <- ifelse(group == groupA, "red", "blue")

png("PCA_T2D_vs_ND_colored.png", width = 1200, height = 900, res = 130)
plot(pca$x[,1], pca$x[,2],
     pch = 19, col = cols,
     xlab = paste0("PC1 (", round(100*summary(pca)$importance[2,1],1), "%)"),
     ylab = paste0("PC2 (", round(100*summary(pca)$importance[2,2],1), "%)"),
     main = "PCA: T2D vs ND")
legend("topright", legend = levels(group),
       col = c("blue", "red"), pch = 19, bty = "n")
dev.off()
cat("✅ Saved: PCA_T2D_vs_ND_colored.png\n")

# 8) ROC curves for top markers (optional; requires pROC from CRAN)
roc_ok <- TRUE
if (!requireNamespace("pROC", quietly = TRUE)) {
  try(install.packages("pROC"), silent = TRUE)
}
if (!requireNamespace("pROC", quietly = TRUE)) {
  roc_ok <- FALSE
  cat("⚠️ pROC not available; skipping ROC curves.\n")
}
if (roc_ok) {
  suppressPackageStartupMessages(library(pROC))
  # pick top 3 significant genes by padj
  top3 <- head(res_df$gene[res_df$padj < 0.05], 3)
  if (length(top3) == 0) top3 <- head(res_df$gene, 3)
  png("ROC_top_genes_T2D_vs_ND.png", width = 1200, height = 900, res = 130)
  plot.new()
  plot.window(xlim = c(0,1), ylim = c(0,1))
  axis(1); axis(2); box()
  title(main = "ROC Curves: Top Genes", xlab = "1 - Specificity", ylab = "Sensitivity")
  # plot each ROC
  for (i in seq_along(top3)) {
    gene <- top3[i]
    roc_obj <- roc(response = group, predictor = as.numeric(log_counts[gene, ]), quiet = TRUE)
    lines(1 - roc_obj$specificities, roc_obj$sensitivities, lwd = 2)
    text(0.6, 0.2 - 0.08*(i-1), labels = paste0(gene, " AUC=", round(auc(roc_obj),3)))
  }
  # diagonal
  segments(0,0,1,1,lty=2)
  dev.off()
  cat("✅ Saved: ROC_top_genes_T2D_vs_ND.png\n")
}

# 9) Export a concise shortlist of biomarkers
short <- subset(res_df, padj < 0.05 & abs(log2FC) > 1)
write.csv(short, "Biomarker_shortlist_T2D_vs_ND.csv", row.names = FALSE)
cat("✅ Saved: Biomarker_shortlist_T2D_vs_ND.csv (padj<0.05 & |log2FC|>1)\n")



