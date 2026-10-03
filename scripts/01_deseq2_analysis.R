suppressPackageStartupMessages({
    library(DESeq2)
    library(ggplot2)
})

message("[1/5] Loading coldata metadata...")
coldata <- read.delim("data/metadata/coldata.tsv", stringsAsFactors = FALSE)
rownames(coldata) <- coldata$sample_id
coldata$condition <- factor(coldata$condition, levels = c("Primary", "Metastatic"))

message("[2/5] Loading raw counts matrix...")
counts_raw <- read.delim("data/raw_counts/GSE278174_rawcounts.txt.gz", stringsAsFactors = FALSE)

# ترکیب Gene Symbol و Ensembl ID برای اطمینان از یکتایی و خوانایی شناسه ژن‌ها
gene_ids <- paste0(counts_raw$Gene, " (", counts_raw$Ensembl_ID, ")")
count_matrix <- as.matrix(counts_raw[, -(1:2)])
rownames(count_matrix) <- gene_ids

# همگام‌سازی و چینش دقیق ستون‌ها با سطر‌های متادیتا
count_matrix <- count_matrix[, coldata$sample_id]
storage.mode(count_matrix) <- "integer"

stopifnot(all(colnames(count_matrix) == rownames(coldata)))
message(sprintf("Loaded count matrix: %d genes across %d samples.", nrow(count_matrix), ncol(count_matrix)))

message("[3/5] Constructing DESeq2 dataset and filtering low counts...")
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = coldata,
                              design = ~ condition)

keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
message(sprintf("Retained %d genes after filtering low counts.", nrow(dds)))

message("[4/5] Running DESeq2 differential expression...")
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Metastatic", "Primary"))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$pvalue), ]

write.csv(res_df, "results/deseq2/differential_expression_Metastatic_vs_Primary.csv", row.names = FALSE)
message("DEG results saved to results/deseq2/differential_expression_Metastatic_vs_Primary.csv")

# خلاصه آماری ژن‌های معنادار
sig_up <- sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE)
sig_down <- sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE)
message(sprintf("Significance (padj < 0.05 & |log2FC| > 1): Up-regulated in Metastatic = %d | Down-regulated = %d", sig_up, sig_down))

message("[5/5] Generating QC & Visualization plots (PCA & Volcano)...")
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = c("condition", "tissue"), returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
    geom_point(size = 3.5, alpha = 0.85) +
    scale_color_manual(values = c("Primary" = "#1f77b4", "Metastatic" = "#d62728")) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle("PCA: GSE278174 Clear Cell RCC (Primary vs Metastatic)") +
    theme_bw(base_size = 13) +
    theme(legend.position = "right", plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("results/plots/pca_plot_primary_vs_metastatic.png", plot = pca_plot, width = 7.5, height = 5.5, dpi = 300)

# Volcano Plot
res_df$Significance <- "Not Significant"
res_df$Significance[res_df$padj < 0.05 & res_df$log2FoldChange > 1] <- "Up in Metastatic"
res_df$Significance[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down in Metastatic"
res_df$Significance <- factor(res_df$Significance, levels = c("Not Significant", "Down in Metastatic", "Up in Metastatic"))

volcano_plot <- ggplot(res_df[!is.na(res_df$padj), ], aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
    geom_point(alpha = 0.5, size = 1.6) +
    scale_color_manual(values = c("Not Significant" = "grey75", "Down in Metastatic" = "#1f77b4", "Up in Metastatic" = "#d62728")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.5) +
    labs(title = "Volcano: Metastatic vs Primary Clear Cell RCC",
         x = expression(log[2] ~ "Fold Change"),
         y = expression(-log[10] ~ "Adjusted P-value")) +
    theme_bw(base_size = 13) +
    theme(legend.position = "top", plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("results/plots/volcano_plot.png", plot = volcano_plot, width = 7, height = 5.5, dpi = 300)
message("Plots successfully saved in results/plots/ directory.")
