#!/usr/bin/env Rscript

suppressWarnings({
  
  suppressPackageStartupMessages({
    library(optparse)
    library(dplyr)
    library(bambu)
  })
  options(dplyr.summarise.inform = FALSE)
  
  option_list = list(
    make_option(c("-b", "--bam"), type="character", default=NULL,
                help="merged bam file"),
    make_option(c("-f", "--fasta"), type="character", default=NULL,
                help="genome reference fasta"),
    make_option(c("-t", "--annotation"), type="character", default=NULL,
                help="reference annotations gtf"),
    make_option(c("-n", "--ndr"), type="numeric", default=NULL,
                help="NDR"),
    make_option(c("-g", "--genefraction"), type="numeric", default=NULL,
                help="minimum readFractionByGene"),
    make_option(c("-o", "--output_dir"), type="character", default=NULL,
                help="path to output directory")
  )
  
  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)
  
  if (is.null(opt$bam)){
    print_help(opt_parser)
    stop("At least one argument must be supplied (input file).n", call.=FALSE)
  }
  
  test.bam <- (opt$bam)
  fa.file <- (opt$fasta)
  gtf.file <- (opt$annotation)
  
  ndr_param <- (opt$ndr)
  genefraction_param <- (opt$genefraction)
  
  bambuAnnotations <- prepareAnnotations(gtf.file)
  
  se <- bambu(reads = test.bam, annotations = bambuAnnotations, genome = fa.file, NDR=ndr_param, opt.discovery=list(min.readFractionByGene=genefraction_param))
  
  # output bambu isoforms
  writeBambuOutput(se, opt$output_dir)
  
  # output tx classes
  tx_data <- as.data.frame(mcols(se))
  tx_data <- apply(tx_data,2,as.character)
  write.table(tx_data, paste0(opt$output_dir, "/bambu_tx_classes.txt"), row.names=F, quote=F, sep="\t")
  
})