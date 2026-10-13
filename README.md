# ccRCC GSE278174: Metastatic vs Primary Transcriptomic Analysis

A reproducible bulk RNA-seq differential expression and functional pathway enrichment workflow analyzing clear cell Renal Cell Carcinoma (ccRCC) patient cohorts from NCBI GEO dataset **GSE278174**.

---

## Overview

- **Dataset:** NCBI GEO [GSE278174](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE278174)
- **Organism:** *Homo sapiens*
- **Design:** Metastatic clear cell Renal Cell Carcinoma (ccRCC) vs Primary ccRCC
- **Cohort Size:** 42 samples (19 Primary, 23 Metastatic)

---

## Key Biological Findings

### 1. Differential Expression Analysis (DESeq2)
- **Total Genes Evaluated:** 14,971
- **Significantly Regulated Genes:** 3,460 genes (1,713 Up, 1,747 Down)

### 2. Functional Pathway Enrichment (GO Biological Process)
Top findings include:
- **Nuclear division** (GO:0000280)
- **Organelle fission** (GO:0048285)
- **Mitotic nuclear division** (GO:0140014)

---

## Repository Structure


```text
.
├── data/raw_counts/GSE278174_rawcounts.txt.gz
├── scripts/01_deseq2_analysis.R
├── scripts/02_pathway_enrichment.R
└── results/

```

---

## Pipeline Execution

1. `conda activate r-ccrcc`
2. `Rscript scripts/01_deseq2_analysis.R`
3. `Rscript scripts/02_pathway_enrichment.R`