# Install R dependencies for the CXCL12 regional Hi-C figure script.

cran_pkgs <- c("strawr", "ggplot2", "patchwork", "ragg")
bioc_pkgs <- c("rtracklayer", "GenomicRanges", "IRanges")

missing_cran <- cran_pkgs[!vapply(cran_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran) > 0) {
  install.packages(missing_cran, repos = "https://cloud.r-project.org")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}

missing_bioc <- bioc_pkgs[!vapply(bioc_pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc) > 0) {
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
}

message("Dependency check complete.")
