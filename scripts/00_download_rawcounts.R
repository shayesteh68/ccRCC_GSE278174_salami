#!/usr/bin/env Rscript

url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE278nnn/GSE278174/suppl/GSE278174_rawcounts.txt.gz"
dest <- "data/raw_counts/GSE278174_rawcounts.txt.gz"

dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)

validate_gzip_table <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  header_line <- readLines(con, n = 1L, warn = FALSE)
  if (length(header_line) != 1L || !nzchar(header_line)) return(FALSE)
  header <- strsplit(header_line, "\t", fixed = TRUE)[[1]]
  if (length(header) < 3L || !identical(header[1:2], c("Ensembl_ID", "Gene")) ||
      anyDuplicated(header)) return(FALSE)
  rows <- 0L
  repeat {
    chunk <- readLines(con, n = 10000L, warn = FALSE)
    if (!length(chunk)) break
    rows <- rows + length(chunk)
  }
  rows > 0L
}

download_rawcounts <- function(url, dest) {
  if (file.exists(dest) && file.info(dest)$size > 0) {
    if (!validate_gzip_table(dest)) stop("Existing count file is not a complete gzip count table: ", dest)
    message("Already present; full gzip and table validation passed: ", dest)
    return(invisible(dest))
  }

  tmp <- paste0(dest, ".part")
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  download.file(url, tmp, mode = "wb", method = "libcurl", quiet = FALSE)
  if (!validate_gzip_table(tmp)) stop("Downloaded file is not a complete gzip count table.")
  if (file.exists(dest)) unlink(dest)
  if (!file.rename(tmp, dest)) stop("Could not move verified download to ", dest)
  message("Saved: ", dest)
  invisible(dest)
}

download_rawcounts(url, dest)
