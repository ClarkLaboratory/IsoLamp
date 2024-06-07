#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(optparse) 
})

# importing
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="input JWR output file"),
  make_option(c("-j", "--JAQ"), type="numeric", default=NULL,
              help="minimum JAQ"),
  make_option(c("-o", "--output_prefix"), type="character", default=NULL,
              help="output name")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

reads <- read.csv(opt$input)
min_jaq <- (opt$JAQ)
output_prefix <- (opt$output_prefix)

reads_filt <- reads %>% 
  group_by(id) %>% 
  filter(min(JAQ) > min_jaq) %>% 
  slice_head(n=1)

write.table(reads_filt$id, paste0(output_prefix, "/temp_files/", "reads_above_JAQ_minimum.txt"), col.names=F, row.names=F, quote=F, sep="\t")
#write.table(reads_filt$id, "~/Documents/sirv_benchmarking_JAQ/sirv5_i/temp_files/sirv5_i_reads_above_JAQ_minimum_0.9.txt", col.names=F, row.names=F, quote=F, sep="\t")
