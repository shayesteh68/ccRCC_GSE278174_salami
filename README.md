# End-to-End Bulk RNA-Seq Pipeline: Clear Cell Renal Cell Carcinoma (ccRCC)

[![Platform: WSL/Ubuntu](https://img.shields.io/badge/Platform-WSL%20%2F%20Ubuntu-orange.svg)]()
[![R: DESeq2](https://img.shields.io/badge/R-DESeq2%20v1.42-blue.svg)]()
[![Analysis: Bulk RNA-Seq](https://img.shields.io/badge/Workflow-DEG%20%26%20Enrichment-success.svg)]()
[![Author: Narges Shayesteh](https://img.shields.io/badge/Author-Narges%20Shayesteh-brightgreen.svg)](https://www.linkedin.com/in/narges-shayesteh)

An end-to-end, publication-grade bulk RNA-seq differential expression and functional enrichment workflow analyzing **Metastatic vs. Primary Clear Cell Renal Cell Carcinoma (ccRCC)** samples from NCBI GEO dataset **GSE278174**.

---

## 🔬 Biological Context & Findings

Clear cell Renal Cell Carcinoma (ccRCC) is characterized by prominent angiogenic activity and metastatic propensity. Using transcriptomic profiling of 42 patient samples (**19 Primary ccRCC** vs. **23 Metastatic ccRCC**), this workflow models metastatic progression:

- **Total Quantified Genes:** 14,971
- **Metastatic Up-regulated:** 1,713 genes ($padj < 0.05$, $\log_2\text{FC} > 1.0$)
- **Metastatic Down-regulated:** 1,747 genes ($padj < 0.05$, $\log_2\text{FC} < -1.0$)
- **Key Pathway Drivers:** Hallmark Hypergeometric Over-Representation Analysis identified significant enrichment of the **Epithelial-Mesenchymal Transition (EMT)** axis ($p = 0.0087$), along with **Angiogenesis** and **Interferon Gamma/Inflammatory Signaling**, highlighting stromal remodeling and immune shifts during distant metastasis.

---

## 📊 Publication-Grade Visualizations

| Principal Component Analysis (PCA) | Transcriptomic Volcano Plot | Hallmark Pathway Enrichment |
| :---: | :---: | :---: |
| ![PCA Plot](results/plots/pca_plot_primary_vs_metastatic.png) | ![Volcano Plot](results/plots/volcano_plot.png) | ![Pathways](results/plots/hallmark_pathways_up.png) |

---

## 📂 Project Structure
```text
ccRCC_GSE278174_salami/
├── data/
│   ├── metadata/
│   │   ├── full_metadata.tsv       # NCBI GEO GSM and phenotype mapping
│   │   └── coldata.tsv             # Curated design table (condition: Primary vs Metastatic)
│   └── rawcounts/
│       └── GSE278174_rawcounts.tsv # Raw expression count matrix (14,971 genes x 42 samples)
├── results/
│   ├── deseq2/
│   │   ├── differential_expression_Metastatic_vs_Primary.csv # Complete DEG statistics
│   │   └── pathway_hallmark_metastatic_up.csv                # Pathway enrichment scores
│   └── plots/
│       ├── pca_plot_primary_vs_metastatic.png              # Variance-stabilized sample clustering (300 DPI)
│       ├── volcano_plot.png          # Highlighting top dysregulated genes (300 DPI)
│       └── hallmark_pathways_up.png  # Top MSigDB Hallmark terms (300 DPI)
├── scripts/
│   ├── 01_deseq2_analysis.R          # DESeq2 modeling, VST normalization, and DEG plotting
│   └── 02_pathway_enrichment.R       # Offline MSigDB Hallmark Fisher exact enrichment
├── .gitignore
└── README.md
⚙️ Environment Setup & Reproducibility
This workflow was executed in an isolated Conda R environment configured for reproducible biostatistical computation.

bash
# Clone the repository
git clone https://github.com/shayesteh68/ccRCC_GSE278174_salami.git
cd ccRCC_GSE278174_salami

# Activate the dedicated R environment
conda activate r_deseq_env

# Step 1: Run DESeq2 differential expression modeling & QC plots
Rscript scripts/01_deseq2_analysis.R

# Step 2: Run Hallmark pathway enrichment analysis
Rscript scripts/02_pathway_enrichment.R
👩‍🔬 Author & Contact
Narges Shayesteh

Role: Bioinformatician & RNA-Seq Data Specialist
LinkedIn: linkedin.com/in/narges-shayesteh
GitHub: @shayesteh68
