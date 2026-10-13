#!/usr/bin/env Rscript
# ==============================================================================
# Functional Pathway Enrichment Analysis (GO-BP & KEGG)
# Dataset: GSE278174 (Clear Cell Renal Cell Carcinoma: Metastatic vs Primary)
# Pipeline: clusterProfiler + org.Hs.eg.db
# ==============================================================================

suppressPackageStartupMessages({
    library(clusterProfiler)
    library(org.Hs.eg.db)
    library(enrichplot)
    library(ggplot2)
})

deg_file <- "results/deseq2/differential_expression_Metastatic_vs_Primary.csv"
if (!file.exists(deg_file)) {
    stop("Input DEG file not found: ", deg_file)
}

message("[1/5] Loading differential expression results...")
deg_df <- read.csv(deg_file, stringsAsFactors = FALSE)

# Clean and filter significant DEGs
deg_df <- deg_df[!is.na(deg_df$padj) & !is.na(deg_df$log2FoldChange), ]
deg_up <- deg_df[deg_df$padj < 0.05 & deg_df$log2FoldChange > 1, ]
deg_down <- deg_df[deg_df$padj < 0.05 & deg_df$log2FoldChange < -1, ]

message(sprintf("Identified %d Up-regulated and %d Down-regulated genes.", nrow(deg_up), nrow(deg_down)))

message("[2/5] Mapping Ensembl IDs to Entrez IDs and Gene Symbols...")
universe_ids <- unique(deg_df$gene)
gene_mapping <- bitr(
    deg_up$gene,
    fromType = "ENSEMBL",
    toType = c("ENTREZID", "SYMBOL"),
    OrgDb = org.Hs.eg.db
)

mapped_entrez_up <- unique(gene_mapping$ENTREZID[!is.na(gene_mapping$ENTREZID)])
message(sprintf("Successfully mapped %d unique Entrez IDs for Up-regulated genes.", length(mapped_entrez_up)))

universe_mapping <- bitr(
    universe_ids,
    fromType = "ENSEMBL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
)
universe_entrez <- unique(universe_mapping$ENTREZID[!is.na(universe_mapping$ENTREZID)])

message("[3/5] Performing Gene Ontology (Biological Process) enrichment...")
ego_bp <- enrichGO(
    gene = mapped_entrez_up,
    universe = universe_entrez,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    readable = TRUE
)

if (!is.null(ego_bp) && nrow(as.data.frame(ego_bp)) > 0) {
    write.csv(as.data.frame(ego_bp), "results/pathway/go_bp_enrichment_metastatic_up.csv", row.names = FALSE)
    p_go <- dotplot(ego_bp, showCategory = 15, title = "GO: Biological Processes (Up in Metastatic ccRCC)") +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
    ggsave("results/plots/go_bp_enrichment_dotplot.png", plot = p_go, width = 8.5, height = 6.5, dpi = 300)
    message("GO BP enrichment completed and saved.")
} else {
    message("No significant GO BP terms found at the selected threshold.")
}

message("[4/5] Performing KEGG Pathway enrichment...")
# Using tryCatch to handle potential network connectivity issues with the KEGG API
ekegg <- tryCatch({
    enrichKEGG(
        gene = mapped_entrez_up,
        universe = universe_entrez,
        organism = "hsa",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
    )
}, error = function(e) {
    message("KEGG query failed: ", e$message)
    return(NULL)
})

if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    ekegg <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    write.csv(as.data.frame(ekegg), "results/pathway/kegg_enrichment_metastatic_up.csv", row.names = FALSE)
    p_kegg <- dotplot(ekegg, showCategory = 15, title = "KEGG Pathways (Up in Metastatic ccRCC)") +
        theme_bw(base_size = 12) +
        theme(plot.title = element_text(face = "bold", hjust = 0.5))
    ggsave("results/plots/kegg_enrichment_dotplot.png", plot = p_kegg, width = 8.5, height = 6.5, dpi = 300)
    message("KEGG pathway enrichment completed and saved.")
} else {
    message("No significant KEGG pathways found or query could not be completed.")
}

message("[5/5] Pathway enrichment workflow finished successfully.")
