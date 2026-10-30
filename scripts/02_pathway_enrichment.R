#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

deg_file <- "results/deseq2/differential_expression_Metastatic_vs_Primary.csv"
if (!file.exists(deg_file)) stop("Input results not found: ", deg_file)
dir.create("results/pathway", recursive = TRUE, showWarnings = FALSE)
dir.create("results/plots", recursive = TRUE, showWarnings = FALSE)

deg <- read.csv(deg_file, stringsAsFactors = FALSE)
deg <- deg[!is.na(deg$gene) & !is.na(deg$log2FoldChange), ]
if (!nrow(deg)) stop("No usable genes in differential-expression results.")

all_map <- bitr(unique(deg$gene), fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
universe <- unique(all_map$ENTREZID)
if (!length(universe)) stop("No Ensembl IDs mapped to Entrez IDs.")

gene_sets <- list(
  metastatic_up = unique(deg$gene[deg$padj < 0.05 & deg$log2FoldChange > 1]),
  metastatic_down = unique(deg$gene[deg$padj < 0.05 & deg$log2FoldChange < -1])
)

for (direction in names(gene_sets)) {
  if (!length(gene_sets[[direction]])) {
    message(direction, ": no significant genes; skipping enrichment.")
    next
  }
  mapped <- bitr(gene_sets[[direction]], fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  genes <- unique(mapped$ENTREZID)
  if (length(genes) < 2L) {
    message(direction, ": fewer than two mapped genes; skipping enrichment.")
    next
  }

  go <- enrichGO(
    gene = genes,
    universe = universe,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
  )
  if (!is.null(go) && nrow(as.data.frame(go))) {
    write.csv(as.data.frame(go), file.path("results/pathway", paste0("go_bp_", direction, ".csv")), row.names = FALSE)
    ggsave(file.path("results/plots", paste0("go_bp_", direction, ".png")),
      dotplot(go, showCategory = 15, title = paste("GO BP:", direction)), width = 8, height = 6, dpi = 300)
  } else {
    message("GO BP ", direction, ": no terms passed the configured significance cutoffs.")
  }

  kegg_failed <- FALSE
  kegg <- tryCatch(
    enrichKEGG(gene = genes, universe = universe, organism = "hsa", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2),
    error = function(e) {
      kegg_failed <<- TRUE
      message("KEGG ", direction, " failed: ", conditionMessage(e))
      NULL
    }
  )
  if (!is.null(kegg) && nrow(as.data.frame(kegg))) {
    kegg <- setReadable(kegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    write.csv(as.data.frame(kegg), file.path("results/pathway", paste0("kegg_", direction, ".csv")), row.names = FALSE)
    ggsave(file.path("results/plots", paste0("kegg_", direction, ".png")),
      dotplot(kegg, showCategory = 15, title = paste("KEGG:", direction)), width = 8, height = 6, dpi = 300)
  } else if (!kegg_failed) {
    message("KEGG ", direction, ": no terms passed the configured significance cutoffs.")
  }
}

message("Direction-specific paired pathway enrichment completed.")
