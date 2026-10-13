# ccRCC GSE278174: Metastatic vs Primary Transcriptomic Analysis

A reproducible bulk RNA-seq differential expression and functional enrichment workflow analyzing clear cell Renal Cell Carcinoma (ccRCC) patient cohorts from NCBI GEO dataset **GSE278174**.

---

## Overview

- **Dataset:** NCBI GEO GSE278174
- **Organism:** *Homo sapiens*
- **Design:** Comparison between Metastatic ccRCC tissue and Primary ccRCC tissue
- **Cohort Size:** 42 samples (19 Primary, 23 Metastatic)
- **Quantification:** Gene-level raw read counts (`GSE278174_rawcounts.txt.gz`) across 14,971 genes

---

## Key Biological Findings

### 1. Differential Expression Analysis (DESeq2)
- Total genes analyzed: **14,971**
- Significantly altered genes (|log2FC| >= 1, padj < 0.05): **3,460 genes**
  - **Upregulated in Metastatic:** 1,713 genes
  - **Downregulated in Metastatic:** 1,747 genes

### 2. Functional Pathway Enrichment (GO Biological Process)
Over-representation analysis of upregulated metastatic genes reveals a pronounced activation of cell cycle and mitotic machinery:
- **Nuclear division** (GO:0000280) — padj = 6.56e-5 (Fold Enrichment: 1.89)
- **Organelle fission** (GO:0048285) — padj = 1.65e-4 (Fold Enrichment: 1.79)
- **Mitotic nuclear division** (GO:0140014) — padj = 3.14e-4 (Fold Enrichment: 1.99)
- **Chromosome segregation** (GO:0007059) — padj = 3.29e-4 (Fold Enrichment: 1.79)
- **Regulation of chromosome separation** (GO:1905818) — padj = 1.73e-3 (Fold Enrichment: 2.84)

---

## Pipeline Execution

### Step 1: Environment Setup

```bash
conda activate r-ccrcc

```

### Step 2: Differential Expression Analysis

```bash
Rscript scripts/01_deseq2_analysis.R

```

### Step 3: Functional Enrichment Analysis

```bash
Rscript scripts/02_pathway_enrichment.R

```
