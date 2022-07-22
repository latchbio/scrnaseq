suppressPackageStartupMessages({
  library(alevinQC)
})

args <- commandArgs(trailingOnly = TRUE)

map_dir <- file.path(args[1])
permit_dir <- file.path(args[2])
quant_dir <- file.path(args[3])
sample_name <- args[4]
report_name <- paste(sample_name, "_alevinQC.html", sep = "")

alevinFryQCReport(
  mapDir = map_dir,
  permitDir = permit_dir,
  quantDir = quant_dir,
  sampleId = sample_name,
  outputFile = report_name,
  outputFormat = "html_document"
)
