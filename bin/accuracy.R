#!/usr/bin/env Rscript

# Function for importing BAM file and calc accuracy 
import_bam_get_accuracy <- function(bamfile) {
  bam <- GenomicAlignments::readGAlignments(bamfile, 
                                            use.names = TRUE,
                                            param = ScanBamParam(tag = c("NM", "AS", "tp"),
                                                                 what = c("qname","flag","mapq")))
  
  message("Imported bam file")
  
  # Summarise CIGAR strings
  cigar_table <- cigarOpTable(bam@cigar)
  
  # CIGAR types 
  col_names_extract <- c("M", "X", "=", "I", "D", "S", "H")
  
  # Add summarised CIGAR strings
  mcols(bam)[paste0("nbr", col_names_extract)] <- mapply(function(col) cigar_table[, col], col_names_extract)
  
  message("Summarised CIGAR strings")
  
  bam_data <- as.data.table(bam %>% setNames(NULL), stringsAsFactors = FALSE) %>%
    dplyr::select(-cigar, -njunc)
  
  # eqx flag
  bam_data <- bam_data %>% 
    dplyr::mutate(read_accuracy=(nbrX+nbr.+nbrI+nbrD-NM)/(nbrX+nbr.+nbrI+nbrD))
  
  # no eqx flag
  #bam_data <- bam_data %>% 
    #dplyr::mutate(read_accuracy=(nbrM+nbrI+nbrD-NM)/(nbrM+nbrI+nbrD))
  
  return(bam_data)
  
}


# Function for plotting accuracy
plot_accuracy <- function(bam_primary, output) {
  pdf(paste0(output, "_read_accuracy.pdf"), width=6, height=6)
  plot <- ggplot(data=bam_primary, aes(x=read_accuracy, y=after_stat(scaled))) +
    geom_density(alpha = 0.4, show.legend = FALSE, fill="steelblue3") +
    theme_bw(base_size=14) +
    xlim(0.8,1) +
    xlab("Read accuracy") +
    ylab("Density")
  suppressMessages(print(plot))
  dev.off()
}

# main function
main <- function() {
  
  # inputs from command line
  args <- commandArgs(trailingOnly = TRUE)
  bamfile <- args[1]
  output <- args[2]
  
  suppressPackageStartupMessages({
    library(GenomicAlignments)
    library(dplyr, warn.conflicts = FALSE)
    library(tibble)
    library(tidyr)
    library(ggplot2)
    library(viridis)
    library(data.table)
  })
  
  options(dplyr.summarise.inform = FALSE)
  
  # call functions
  bam_tidy <- import_bam_get_accuracy(bamfile)
  
  # plotting functions
  plot_accuracy(bam_tidy, output)
  
  message("Complete")
  
}

suppressWarnings(
  main())