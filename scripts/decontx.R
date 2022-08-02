suppressPackageStartupMessages({
  library(zellkonverter)
  library(celda)
})

args <- commandArgs(trailingOnly = TRUE)

h5ad <- file.path(args[1])
sce <- readH5AD(h5ad)
sce <- decontX(sce, assayName = "X")

# Extraneous data 3xes the size of the h5ad
reducedDims(sce)$decontX_UMAP <- NULL

assays(sce)$X <- assays(sce)$decontXcounts
assays(sce)$decontXcounts <- NULL

metadata(sce)$decontX <- NULL
writeH5AD(sce, file = args[1])
