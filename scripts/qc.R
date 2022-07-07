suppressPackageStartupMessages({
  library(alevinQC)
})


args = commandArgs(trailingOnly=TRUE)

mapD = file.path(args[0])
permitD = file.path(args[1])
quantD = file.path(args[2])
name = args[3]
reportName = "alevinQC.html"

alevinFryQCReport(
  mapDir = mapD,
	permitDir = permitD,
  quantDir = quantD,
	sampleId = name,
	outputFile = reportName,
	outputFormat="html_document")
