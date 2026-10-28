#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

metadata_path <- "data/metadata/coldata.tsv"
patient_map_path <- "data/metadata/patient_map.tsv"
counts_path <- "data/raw_counts/GSE278174_rawcounts.txt.gz"
out_dir <- "results/deseq2"
plot_dir <- "results/plots"

for (path in c(metadata_path, patient_map_path, counts_path)) {
  if (!file.exists(path)) stop("Required input is missing: ", path)
}
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

coldata <- read.delim(metadata_path, stringsAsFactors = FALSE, check.names = FALSE)
patient_map <- read.delim(patient_map_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("sample_id", "condition") %in% names(coldata))) {
  stop("coldata.tsv must contain sample_id and condition.")
}
if (!all(c("sample_id", "patient_id") %in% names(patient_map))) {
  stop("patient_map.tsv must contain sample_id and patient_id.")
}
if (anyDuplicated(coldata$sample_id) || anyDuplicated(patient_map$sample_id)) {
  stop("Sample IDs must be unique in both metadata files.")
}
if (anyNA(patient_map$patient_id) || any(!nzchar(trimws(patient_map$patient_id)))) {
  stop("Fill in a study-confirmed patient_id for every sample in patient_map.tsv.")
}
if (!setequal(coldata$sample_id, patient_map$sample_id)) {
  stop("patient_map.tsv must contain exactly the sample IDs in coldata.tsv.")
}
coldata$patient_id <- patient_map$patient_id[match(coldata$sample_id, patient_map$sample_id)]
coldata$condition <- factor(coldata$condition, levels = c("Primary", "Metastatic"))
if (anyNA(coldata$condition)) stop("condition values must be Primary or Metastatic.")

pair_table <- table(coldata$patient_id, coldata$condition)
paired_patients <- rownames(pair_table)[
  pair_table[, "Primary"] > 0 & pair_table[, "Metastatic"] > 0
]
if (length(paired_patients) < 2L) {
  stop("At least two patients with both conditions are required for paired DE.")
}
coldata <- coldata[coldata$patient_id %in% paired_patients, , drop = FALSE]
coldata$patient_id <- droplevels(factor(coldata$patient_id))
rownames(coldata) <- coldata$sample_id

counts_raw <- read.delim(counts_path, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("Ensembl_ID", "Gene") %in% names(counts_raw))) {
  stop("Count matrix must contain Ensembl_ID and Gene columns.")
}
missing_samples <- setdiff(coldata$sample_id, names(counts_raw))
if (length(missing_samples)) {
  stop("Samples missing from counts: ", paste(missing_samples, collapse = ", "))
}
count_matrix <- as.matrix(counts_raw[, coldata$sample_id, drop = FALSE])
storage.mode(count_matrix) <- "numeric"
if (anyNA(count_matrix) || any(!is.finite(count_matrix)) ||
    any(count_matrix < 0) || any(count_matrix != floor(count_matrix))) {
  stop("Raw count columns must contain finite, non-negative integers.")
}
storage.mode(count_matrix) <- "integer"
gene_ids <- sub("\\..*$", "", counts_raw$Ensembl_ID)
fallback <- is.na(gene_ids) | !nzchar(gene_ids)
gene_ids[fallback] <- counts_raw$Gene[fallback]
rownames(count_matrix) <- make.unique(gene_ids)
if (!identical(colnames(count_matrix), rownames(coldata))) {
  stop("Count columns and metadata rows are not in the same order.")
}

input_tumor_samples <- nrow(coldata)
sample_key <- paste(coldata$patient_id, coldata$condition, sep = "__")
sample_levels <- unique(sample_key)
aggregated_counts <- vapply(
  sample_levels,
  function(key) rowSums(count_matrix[, sample_key == key, drop = FALSE]),
  numeric(nrow(count_matrix))
)
colnames(aggregated_counts) <- sample_levels
if (anyNA(aggregated_counts) || any(!is.finite(aggregated_counts)) ||
    any(aggregated_counts < 0) || any(aggregated_counts > .Machine$integer.max) ||
    any(aggregated_counts != floor(aggregated_counts))) {
  stop("Patient-condition aggregate counts must be non-negative integers within R's integer range.")
}
storage.mode(aggregated_counts) <- "integer"

first_sample <- match(sample_levels, sample_key)
coldata <- data.frame(
  sample_id = sample_levels,
  patient_id = as.character(coldata$patient_id[first_sample]),
  condition = as.character(coldata$condition[first_sample]),
  n_input_tumor_samples = vapply(sample_levels, function(key) sum(sample_key == key), integer(1)),
  stringsAsFactors = FALSE
)
rownames(coldata) <- coldata$sample_id
coldata$patient_id <- factor(coldata$patient_id)
coldata$condition <- relevel(factor(coldata$condition, levels = c("Primary", "Metastatic")), ref = "Primary")
if (!identical(colnames(aggregated_counts), rownames(coldata))) {
  stop("Aggregated count columns and patient-condition metadata are out of order.")
}
design_matrix <- model.matrix(~ patient_id + condition, data = coldata)
if (qr(design_matrix)$rank < ncol(design_matrix)) {
  stop("Paired patient-condition design is not full rank; inspect patient_map.tsv and conditions.")
}

dds <- DESeqDataSetFromMatrix(
  countData = aggregated_counts,
  colData = coldata,
  design = ~ patient_id + condition
)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep, ]
if (nrow(dds) < 2L) stop("Too few genes remain after low-count filtering.")
dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Metastatic", "Primary"))
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df[order(res_df$padj, na.last = TRUE), , drop = FALSE]
write.csv(
  res_df,
  file.path(out_dir, "differential_expression_Metastatic_vs_Primary.csv"),
  row.names = FALSE
)

sig_up <- sum(res$padj < 0.05 & res$log2FoldChange > 1, na.rm = TRUE)
sig_down <- sum(res$padj < 0.05 & res$log2FoldChange < -1, na.rm = TRUE)
writeLines(c(
  "Paired DESeq2 analysis",
  paste("Patients with both conditions:", length(paired_patients)),
  paste("Original tumor samples retained:", input_tumor_samples),
  paste("Patient-condition aggregates:", nrow(coldata)),
  "Multiple metastatic lesions from one patient are summed within that patient-condition.",
  paste("Genes retained:", nrow(dds)),
  paste("Significant genes (padj < 0.05 and |log2FC| > 1):", sig_up + sig_down),
  paste("Up in metastatic:", sig_up),
  paste("Down in metastatic:", sig_down),
  "Design: ~ patient_id + condition",
  "Contrast: Metastatic versus Primary"
), file.path(out_dir, "summary.txt"))

vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = c("condition", "patient_id"), returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))
pca <- ggplot(pca_data, aes(PC1, PC2, color = condition, shape = patient_id)) +
  geom_point(size = 3, alpha = 0.85) +
  labs(
    title = "Paired ccRCC samples: primary and metastatic",
    x = paste0("PC1: ", percent_var[1], "% variance"),
    y = paste0("PC2: ", percent_var[2], "% variance")
  ) +
  theme_bw()
ggsave(file.path(plot_dir, "pca_paired_primary_vs_metastatic.png"), pca, width = 8, height = 6, dpi = 300)

volcano_data <- res_df[!is.na(res_df$padj), , drop = FALSE]
volcano_data$group <- "Not significant"
volcano_data$group[volcano_data$padj < 0.05 & volcano_data$log2FoldChange > 1] <- "Up in metastatic"
volcano_data$group[volcano_data$padj < 0.05 & volcano_data$log2FoldChange < -1] <- "Down in metastatic"
volcano <- ggplot(volcano_data, aes(log2FoldChange, -log10(pmax(padj, .Machine$double.xmin)), color = group)) +
  geom_point(alpha = 0.55, size = 1.3) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(
    title = "Paired differential expression: metastatic versus primary",
    x = "log2 fold change",
    y = "-log10 adjusted p-value"
  ) +
  theme_bw()
ggsave(file.path(plot_dir, "volcano_paired_primary_vs_metastatic.png"), volcano, width = 8, height = 6, dpi = 300)

message("Paired DE completed for ", length(paired_patients), " patients from ", input_tumor_samples, " tumor samples.")
