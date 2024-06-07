#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})

options(dplyr.summarise.inform = FALSE)

option_list = list(
  make_option(c("-e", "--ens_var"), type="character", default=NULL,
              help="ENS_ID variable"),
  make_option(c("-m", "--proportion_min"), type="character", default=NULL,
              help="Proportion threshold"),
  make_option(c("-s", "--sample_min"), type="character", default=NULL,
              help="Sample threshold"),
  make_option(c("-o", "--output_path"), type="character", default=NULL,
              help="output path")
)

## write code to skip folder if there is no quant.sf inside ##


suppressWarnings({

  opt_parser = OptionParser(option_list=option_list)
  opt = parse_args(opt_parser)

  gene_id <- (opt$ens_var)
  proportion_min <- as.numeric(opt$proportion_min)/1000000
  samples_min <- as.numeric(opt$sample_min)
  outdir <- (opt$output_path)
  
  # get sample ids from within directory
  all_folders <- c(list.files(paste0(outdir,"/updated_transcriptome/salmon_quants/")))
  # get sample ids from within directory
  sample_files <- list.files(path = paste0(outdir,"/updated_transcriptome/salmon_quants/"), pattern = "quant.sf", recursive = TRUE, full.names = FALSE)
  
  if (length(sample_files) < length(all_folders)) {
    print("WARNING: not every sample was quantified with salmon.")
    print("If this was unexpected, check directories for missing quant.sf files: outdir/updated_transcriptome/salmon_quants")
  }
  
  import_quants_function <- function(sample_id) {
    # get paths of file to import
    pathtofile <- paste0(outdir,"/updated_transcriptome/salmon_quants/", sample_id)
    
    df <- read.table(pathtofile, header=T)

    # keep only transcript_id and read counts columns
    df <- df[, c("Name", "NumReads")]
    # replace all "." characters with a "_" in case IsoVis breaks
    sample_id <- gsub("\\.","_",sample_id)
    # rename cols to transcript_id and the sample id
    colnames(df) <- c("transcript_id", paste0(sample_id))
    rownames(df) <- df[,1]
    df[,1] <- NULL
    
    return(df)
  }
  
  # import all files and combine dfs
  all_samples <- lapply(sample_files, import_quants_function)
  
  combined_counts <- do.call("cbind", all_samples)
  
  # format column names
  column_names <- c(colnames(combined_counts))
  split_names <- strsplit(column_names, "/")
  sample_names <- sapply(split_names, function(x) x[1])
  
  colnames(combined_counts) <- sample_names
  
  # remove version numbers
  rownames(combined_counts) <- gsub("\\..*","", rownames(combined_counts))

  # get order of columns
  new_col_order <- c("transcript_id", paste0(colnames(combined_counts)))
  # convert row names to column called transcript_id
  combined_counts$transcript_id <- rownames(combined_counts)
  # remove rownames
  rownames(combined_counts) <- NULL
  # order df
  combined_counts <- combined_counts[, new_col_order]
  
  # proportions first pass
  
  # store transcript_id as vector
  txids <- as.vector(combined_counts$transcript_id)
  
  # copy counts to new df
  combined_props <- combined_counts
  combined_props$transcript_id <- NULL
  
  combined_props_names <-c(paste0(colnames(combined_props)))
  # calculate proportions per sample
  combined_props_df <- data.frame(lapply(combined_props, function(x) x / sum(x)))
  colnames(combined_props_df) <- combined_props_names
  # add transcript_id back and order df
  combined_props_df$transcript_id <- txids
  combined_props_df <- combined_props_df[, c("transcript_id", combined_props_names)]
  
  combined_counts <- combined_counts %>% dplyr::arrange(transcript_id)
  combined_props_df <- combined_props_df %>% dplyr::arrange(transcript_id)
 
  
  # apply read count and samples minimum thresholds
  combined_props_df_filt <- combined_props_df %>% 
    dplyr::filter(rowSums(select_if(., is.numeric) >= proportion_min, na.rm = TRUE) >= samples_min)
  
  # proportions second pass
  
  combined_counts_filtered <- combined_counts %>% filter(transcript_id %in% combined_props_df_filt$transcript_id)
  
  # store transcript_id as vector
  txids <- as.vector(combined_counts_filtered$transcript_id)
  
  # copy counts to new df
  combined_props_final <- combined_counts_filtered
  combined_props_final$transcript_id <- NULL
  
  combined_props_names <-c(paste0(colnames(combined_props_final)))
  # calculate proportions per sample
  combined_props_final_df <- data.frame(lapply(combined_props_final, function(x) x / sum(x)))
  colnames(combined_props_final_df) <- combined_props_names
  # add transcript_id back and order df
  combined_props_final_df$transcript_id <- txids
  combined_props_final_df <- combined_props_final_df[, c("transcript_id", combined_props_names)]
  
  combined_counts_filtered <- combined_counts_filtered %>% dplyr::arrange(transcript_id)
  combined_props_final_df <- combined_props_final_df %>% dplyr::arrange(transcript_id)
  
  # calculate total reads from df
  total_reads_remaining <- as.integer(sum(rowSums(select_if(combined_counts_filtered, is.numeric), na.rm = TRUE)))

  # create TPM values
  combined_TPM <- combined_props_final_df %>% mutate_if(is.numeric, ~ . * 1000000)

  # add gene id col
  combined_counts_filtered <- cbind(gene_id = gene_id, combined_counts_filtered)
  combined_props_final_df <- cbind(gene_id = gene_id, combined_props_final_df)
  combined_TPM <- cbind(gene_id = gene_id, combined_TPM)

  # export files
  write.csv(combined_counts_filtered, paste0(outdir, "/", outdir, "_counts.csv"), quote = FALSE, row.names = FALSE)
  write.csv(combined_props_final_df, paste0(outdir, "/", outdir, "_proportions.csv"), quote = FALSE, row.names = FALSE)
  write.csv(combined_TPM, paste0(outdir, "/", outdir, "_TPM_values.csv"), quote = FALSE, row.names = FALSE)
  
  
  write.table(total_reads_remaining, paste0(outdir, "/", "temp_files/", "remaining_read_sum.txt"), quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  
  
})
    