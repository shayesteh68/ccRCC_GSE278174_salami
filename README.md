# ccRCC primary versus metastatic RNA-seq (GSE278174)

This repository contains a bulk RNA-seq analysis workflow for the 42 samples in GEO series GSE278174: 19 primary tumors and 23 metastatic samples from 19 patients. The design is paired by patient, and several patients contributed more than one metastasis.

## Inputs

- `data/metadata/coldata.tsv` — sample IDs, GEO accessions, tissue sites, and conditions.
- `data/metadata/patient_map.tsv` — sample-to-patient mapping from the study's Supplementary Information 3.
- `data/metadata/PATIENT_MAP_PROVENANCE.md` — source fields and validation summary for the mapping.
- `data/raw_counts/GSE278174_rawcounts.txt.gz` — GEO raw counts; retrieve with `Rscript scripts/00_download_rawcounts.R`.

The verified map links every sample to its de-identified patient code using the study's `Library.ID` and `Myriad` fields. The analysis checks that sample IDs match exactly and that the paired design is usable before running differential expression.

## Run

Create the Conda environment, then run from the repository root:

```bash
conda env create -f environment.yml
conda activate ccrcc-rnaseq
Rscript scripts/00_fetch_geo_metadata.R
Rscript scripts/00_download_rawcounts.R
Rscript scripts/01_deseq2_analysis.R
Rscript scripts/02_pathway_enrichment.R
```

The DE model is `~ patient_id + condition`, with metastatic versus primary as the contrast. This implementation is a paired DESeq2 workflow; it is not an exact reproduction of the study's published edgeR/limma analysis. Existing outputs, if present, should be regenerated before citation.

Code in this repository is licensed under the MIT License (see `LICENSE`). Third-party data, software, and study supplements remain subject to their own terms.
