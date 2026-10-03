suppressPackageStartupMessages({
    library(ggplot2)
})

message("[1/3] Reading DEG results...")
res <- read.csv("results/deseq2/differential_expression_Metastatic_vs_Primary.csv", stringsAsFactors = FALSE)

# تفکیک شناسه و Gene Symbol
res$symbol <- sub(" .*", "", res$gene)

# ژن‌های معنادار بیش‌بیان (Up-regulated) و کلیه ژن‌های بک‌گراند
up_genes <- unique(na.omit(res$symbol[res$padj < 0.05 & res$log2FoldChange > 1.0]))
all_genes <- unique(na.omit(res$symbol))

message(sprintf("Identified %d significant up-regulated genes out of %d total background genes.", length(up_genes), length(all_genes)))

# تعریف ژن‌ست‌های کلیدی MSigDB Hallmark بدون نیاز به دانلود خارجی
message("[2/3] Performing Hypergeometric Enrichment Test on MSigDB Hallmark Pathways...")

hallmark_sets <- list(
  EPITHELIAL_MESENCHYMAL_TRANSITION = c(
    "VIM", "FN1", "CDH2", "SNAI1", "SNAI2", "TWIST1", "ZEB1", "ZEB2", "MMP2", "MMP9", 
    "ACTA2", "COL1A1", "COL1A2", "COL3A1", "COL5A1", "COL5A2", "TGFB1", "TGFB2", "TGFB3",
    "SERPINE1", "SPARC", "VCAN", "LAMC2", "COL6A2", "COL6A3", "LOX", "LOXL2"
  ),
  HYPOXIA = c(
    "VEGFA", "SLC2A1", "CA9", "LDHA", "ENO1", "PDK1", "PGK1", "HK2", "BNIP3", "ALDOA",
    "GAPDH", "ADM", "NDRG1", "ANKRD37", "PFKFB3", "EGLN3", "LOX", "CXCR4", "BHLHE40"
  ),
  ANGIOGENESIS = c(
    "VEGFA", "KDR", "FLT1", "ANGPT1", "ANGPT2", "TEK", "PECAM1", "CD34", "ENG", "VWF",
    "CCND1", "PDGFA", "PDGFB", "COL3A1", "POSTN", "LUM", "NRP1", "CXCR4", "SERPINE1"
  ),
  INFLAMMATORY_RESPONSE = c(
    "IL6", "CXCL8", "CCL2", "ICAM1", "VCAM1", "SELE", "CD44", "PTGS2", "IRF1", "IRF7",
    "CXCL1", "CXCL2", "CXCL10", "NFKB1", "RELA", "IL1B", "TNF", "TLR2", "TLR4"
  ),
  TNFA_SIGNALING_VIA_NFKB = c(
    "NFKB1", "NFKBIA", "TNFAIP3", "JUN", "FOS", "JUNB", "IL6", "CXCL1", "CXCL2", "CXCL3",
    "BIRC2", "BIRC3", "RELB", "ICAM1", "CCL2", "CXCL8", "CD83", "DUSP1", "DUSP2", "EGR1"
  ),
  GLYCOLYSIS = c(
    "HK2", "HK1", "GPI", "PFKP", "ALDOA", "ALDOC", "TPI1", "GAPDH", "PGK1", "PGAM1",
    "ENO1", "PKM", "LDHA", "SLC2A1", "G6PD", "PDK1", "PDK3", "PFKFB3", "PGM1"
  ),
  INTERFERON_GAMMA_RESPONSE = c(
    "STAT1", "IRF1", "IRF9", "ISG15", "IFIT1", "IFIT2", "IFIT3", "CXCL9", "CXCL10", 
    "CXCL11", "OAS1", "OAS2", "GBP1", "GBP2", "B2M", "HLA-A", "HLA-B", "CD74"
  ),
  COMPLEMENT_SYSTEM = c(
    "C1R", "C1S", "C3", "C4A", "C4B", "CFB", "CFH", "CFI", "CD55", "CD59", "SERPING1",
    "C1QA", "C1QB", "C1QC", "ITGB2", "ITGAM", "A2M", "CASP1", "CP"
  ),
  COAGULATION = c(
    "F3", "PLAT", "PLAU", "SERPINE1", "SERPINC1", "FGB", "FGG", "PROS1", "PROC", "THBD",
    "VWF", "F10", "F7", "F5", "ANXA1", "ANXA2", "ANXA5", "CD36", "TIMP1"
  ),
  ALLOGRAFT_REJECTION = c(
    "CD3D", "CD3E", "CD3G", "CD4", "CD8A", "CD8B", "IL2RA", "PTPRC", "IFNG", "PRF1",
    "GZMA", "GZMB", "FASLG", "CD28", "ICOS", "CTLA4", "PDCD1", "HLA-DRA", "HLA-DRB1"
  )
)

results_list <- list()

for (name in names(hallmark_sets)) {
  pathway_genes <- hallmark_sets[[name]]
  
  q <- length(intersect(up_genes, pathway_genes))
  m <- length(intersect(all_genes, pathway_genes))
  n <- length(all_genes) - m
  k <- length(up_genes)
  
  p_val <- phyper(q - 1, m, n, k, lower.tail = FALSE)
  
  overlap_genes <- paste(intersect(up_genes, pathway_genes), collapse = ", ")
  
  results_list[[name]] <- data.frame(
    Pathway = gsub("_", " ", name),
    Overlap = q,
    Pathway_Size = m,
    P_value = p_val,
    Genes = overlap_genes,
    stringsAsFactors = FALSE
  )
}

pathway_df <- do.call(rbind, results_list)
pathway_df$padj <- p.adjust(pathway_df$P_value, method = "BH")
pathway_df <- pathway_df[order(pathway_df$P_value), ]

write.csv(pathway_df, "results/deseq2/pathway_hallmark_metastatic_up.csv", row.names = FALSE)
message("Saved pathway results to results/deseq2/pathway_hallmark_metastatic_up.csv")

message("[3/3] Generating publication-grade Pathway Barplot...")
pathway_df$Pathway <- factor(pathway_df$Pathway, levels = rev(pathway_df$Pathway))

p <- ggplot(pathway_df, aes(x = -log10(pmax(padj, 1e-10)), y = Pathway, fill = Overlap)) +
  geom_col(width = 0.7, color = "black") +
  scale_fill_gradient(low = "#ffb3b3", high = "#990000") +
  labs(
    title = "Key Hallmark Pathways Enriched in Metastatic ccRCC",
    subtitle = "GSE278174: Metastatic vs Primary Clear Cell Renal Cell Carcinoma",
    x = expression(-log[10] ~ "(Adjusted P-value)"),
    y = "",
    fill = "Gene Overlap"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5),
    axis.text.y = element_text(face = "bold", color = "black"),
    panel.grid.minor = element_blank()
  )

ggsave("results/plots/hallmark_pathways_up.png", plot = p, width = 8.5, height = 5.2, dpi = 300)
message("Pathway plot successfully saved in results/plots/hallmark_pathways_up.png")

print(pathway_df[, c("Pathway", "Overlap", "Pathway_Size", "P_value", "padj")])
