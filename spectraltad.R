#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
matrix <- args[which(args == "--matrix") + 1]
out_dir <- args[which(args == "--output") + 1]
if (!grepl("/$", out_dir)) {
    out_dir <- paste0(out_dir, "/")
}

library(SpectralTAD)
library(HiCcompare)

cool_mat <- read.table(matrix)
sparse_mats = HiCcompare::cooler2sparse(cool_mat)
spec_tads = lapply(names(sparse_mats), function(x) {
  tryCatch({
    SpectralTAD(sparse_mats[[x]], chr = x, levels = 3)
  }, error = function(e) {
    message(paste("Error (levels=3) with chromosome", x, ":", e$message))
    message(paste("Trying levels=2 for", x))    
    tryCatch({
      SpectralTAD(sparse_mats[[x]], chr = x, levels = 2)
    }, error = function(e2) {
      message(paste("Failed even for levels 1-2 for", x, ":", e2$message))
      return(NULL)
    })
  })
})
for(i in seq_along(spec_tads)) {
  if(!is.null(spec_tads[[i]])) {
    chr_name <- names(sparse_mats)[i]
    for(level in 1:length(spec_tads[[i]])) {
      write.table(spec_tads[[i]][[level]],
                  paste0(out_dir, "SpectralTAD_", chr_name, "_level", level, ".bed"),
                  row.names=FALSE, col.names=TRUE, sep="\t", quote=FALSE)
    }
  }
}
