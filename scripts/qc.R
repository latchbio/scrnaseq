suppressPackageStartupMessages({
  library(alevinQC)
})

args <- commandArgs(trailingOnly = TRUE)

map_dir <- file.path(args[0])
permit_dir <- file.path(args[1])
quant_dir <- file.path(args[2])
name <- args[3]
report_name <- "alevinQC.html"

alevinFryQCReport(
  mapDir = map_dir,
  permitDir = permit_dir,
  quantDir = quant_dir,
  sampleId = name,
  outputFile = report_name,
  outputFormat = "html_document"
)
