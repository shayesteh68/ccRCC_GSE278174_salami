# Clear Cell Renal Cell Carcinoma (ccRCC) Progression & Metastasis Analysis (GSE278174)

## Cohort Summary

This repository contains a complete, reproducible RNA-Seq differential expression analysis of **clear cell renal cell carcinoma (ccRCC)** progression and metastasis.

- **Cohort:** Clinical cohort of Dr. Simpa S. Salami (University of Michigan Rogel Cancer Center)
- **Data source:** NCBI GEO accession [GSE278174](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278174)
- **Samples:** 42 patient samples — **19 Primary ccRCC** vs **23 Metastatic ccRCC**
- **Data type:** RNA-Seq raw counts, **~14,971 genes** profiled
- **Analysis:** DESeq2 differential expression (Metastatic vs Primary)

## Key Results

Differential expression of metastatic vs primary tumors (|log2FC| >= 1, FDR < 0.05):

- **1,713 Up-regulated** genes in metastatic ccRCC
- **1,747 Down-regulated** genes in metastatic ccRCC

Hallmark pathway enrichment was performed **offline** using a Fisher exact test against MSigDB Hallmark gene sets. Top findings:

- **Epithelial-Mesenchymal Transition (EMT):** p = 0.0087
- **Interferon Gamma Response:** p = 0.0381

## Publication Figures

| Sample Segregation (PCA) | Global Transcriptomic Shift (Volcano) | Hallmark Enrichment |
| :---: | :---: | :---: |
| ![PCA](results/plots/pca_plot_primary_vs_metastatic.png) | ![Volcano](results/plots/volcano_plot.png) | ![Hallmark Pathways](results/plots/hallmark_pathways_up.png) |

## Project Structure

```text
ccRCC_GSE278174_salami/
├── data/
│   ├── metadata/
│   │   └── coldata.tsv                  # Sample metadata: 42 samples, condition labels
│   └── rawcounts/
│       └── GSE278174_rawcounts.tsv      # Raw count matrix, ~14,971 genes
├── results/
│   ├── deseq2/
│   │   ├── differential_expression_Metastatic_vs_Primary.csv   # Full DESeq2 results table
│   │   └── pathway_hallmark_metastatic_up.csv                  # Hallmark enrichment results
│   └── plots/
│       ├── pca_plot_primary_vs_metastatic.png                  # PCA of sample segregation
│       ├── volcano_plot.png                                    # Volcano plot of DE genes
│       └── hallmark_pathways_up.png                            # Upregulated hallmark pathways
├── scripts/
│   ├── 00_fetch_geo_metadata.R         # Fetch and verify GEO sample metadata
│   ├── 01_deseq2_analysis.R            # DESeq2 analysis, PCA and volcano plots
│   ├── 02_pathway_enrichment.R         # Offline Hallmark Fisher exact enrichment
│   └── 02_pathway_enrichment.py        # Python variant of enrichment workflow
├── .gitignore
├── logs/
└── README.md
```

## Environment Setup & Reproducibility

```bash
# Clone the repository
git clone https://github.com/shayesteh68/ccRCC_GSE278174_salami.git
cd ccRCC_GSE278174_salami

# Activate the dedicated R/DESeq2 conda environment
conda activate r_deseq_env

# Step 1: Fetch and verify sample metadata (42 samples)
Rscript scripts/00_fetch_geo_metadata.R

# Step 2: DESeq2 differential expression + PCA/Volcano plots
Rscript scripts/01_deseq2_analysis.R

# Step 3: Offline Hallmark pathway enrichment (Fisher exact test)
Rscript scripts/02_pathway_enrichment.R
```

## Author

**Narges Shayesteh**, Bioinformatician & RNA-Seq Data Specialist

- LinkedIn: <https://www.linkedin.com/in/narges-shayesteh>
- GitHub: [@shayesteh68](https://github.com/shayesteh68)
