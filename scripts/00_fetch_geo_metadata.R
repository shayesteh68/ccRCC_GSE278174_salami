suppressPackageStartupMessages({
  library(GEOquery)
})

gse_id <- "GSE278174"
out_rds <- file.path("data", "metadata", paste0(gse_id, "_GSE.rds"))

message("[1/3] Downloading GEO SOFT metadata (lightweight)...")
gse <- getGEO(gse_id, GSEMatrix = TRUE, getGPL = FALSE)

saveRDS(gse, out_rds)
message("[2/3] Saved: ", out_rds)

# Extract pheno if expression matrix exists in GSEMatrix
if (is.list(gse) && length(gse) >= 1) {
  eset <- gse[[1]]
  pheno <- pData(eset)
  write.csv(pheno, file.path("data","metadata", paste0(gse_id, "_pheno.csv")), row.names = FALSE)
  message("[3/3] Wrote phenotype: data/metadata/", gse_id, "_pheno.csv")
} else {
  message("[WARN] GSEMatrix not available as expression set; may be SRA-only. We'll handle in next step.")
}
