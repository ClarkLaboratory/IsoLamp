#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(GenomicAlignments)
  library(dplyr, warn.conflicts = FALSE)
  library(tibble)
  library(tidyr)
  library(data.table)
  library(dplyr)
  library(optparse)
})

# importing
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="input BAM file"),
  make_option(c("-a", "--accuracy"), type="numeric", default=NULL,
              help="minimum read accuracy"),
  make_option(c("-o", "--output_prefix"), type="character", default=NULL,
              help="output name")
)


opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

bamfile <- (opt$input)
min_acc <- (opt$accuracy)
output_prefix <- (opt$output_prefix)

bam <- GenomicAlignments::readGAlignments(bamfile, 
                                            use.names = TRUE,
                                            param = ScanBamParam(tag = c("NM", "AS", "tp"),
                                                                 what = c("qname","flag","mapq")))
# Summarise CIGAR strings
cigar_table <- cigarOpTable(bam@cigar)

# CIGAR types 
col_names_extract <- c("M", "X", "=", "I", "D", "S", "H")
  
# Add summarised CIGAR strings
mcols(bam)[paste0("nbr", col_names_extract)] <- mapply(function(col) cigar_table[, col], col_names_extract)
  
bam_data <- as.data.table(bam %>% setNames(NULL), stringsAsFactors = FALSE) %>%
  dplyr::select(-cigar, -njunc)
  
# eqx flag
bam_data <- bam_data %>% 
  dplyr::mutate(read_accuracy=(nbrX+nbr.+nbrI+nbrD-NM)/(nbrX+nbr.+nbrI+nbrD))

bam_filt <- bam_data[read_accuracy > min_acc]

write.table(bam_filt$qname, paste0(output_prefix, "/temp_files/", "reads_above_accuracy_minimum.txt"), col.names=F, row.names=F, quote=F, sep="\t")






  
  