suppressPackageStartupMessages({
  library(eisaR)
  library(Biostrings)
  library(BSgenome)
  library(stringr)
  library(GenomicFeatures)
})


args = commandArgs(trailingOnly=TRUE)

gtf_path = file.path(args[0])
genome_path = file.path(args[1])
read_length = as.integer(args[2])
flank_trim_length = 5
output_dir = "splici_txome"

make_splici_txome(
  gtf_path=gtf_path,
	genome_path=genome_path,
	read_length=read_length,
	flank_trim_length=flank_trim_length,
	output_dir=output_dir)
