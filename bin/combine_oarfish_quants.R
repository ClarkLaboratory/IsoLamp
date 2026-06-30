#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
})

options(dplyr.summarise.inform = FALSE)

# Command-line options
option_list = list(
  make_option(c("-e", "--ens_var"), type="character", default=NULL,
              help="Gene/ENS_ID variable"),
  make_option(c("-m", "--proportion_min"), type="character", default=NULL,
              help="Minimum proportion threshold (CPM)"),
  make_option(c("-s", "--sample_min"), type="character", default=NULL,
              help="Minimum number of samples for filter"),
  make_option(c("-o", "--output_path"), type="character", default=NULL,
              help="Output directory path")
)

suppressWarnings({

  opt_parser <- OptionParser(option_list = option_list)
  opt <- parse_args(opt_parser)

  gene_id <- opt$ens_var
  proportion_min <- as.numeric(opt$proportion_min)/1e6  # convert CPM to proportion
  samples_min <- as.numeric(opt$sample_min)
  outdir <- opt$output_path

  # Locate all Oarfish .quant files
  quant_dir <- file.path(outdir, "oarfish_quants")
  sample_files <- list.files(path = quant_dir, pattern = "\\.quant$", recursive = FALSE, full.names = TRUE)

  if(length(sample_files) == 0){
    stop("No .quant files found in ", quant_dir)
  }

  # Function to import each .quant file and assign proper sample name
  import_oarfish_quant <- function(file_path){
    df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    df <- df[, c("tname", "num_reads")]  # only keep transcript_id and read counts

    # Extract sample name from file basename
    sample_id <- tools::file_path_sans_ext(basename(file_path))
    sample_id <- gsub("\\.", "_", sample_id)  # replace dots with underscores
    colnames(df) <- c("transcript_id", sample_id)

    rownames(df) <- df$transcript_id
    df$transcript_id <- NULL
    return(df)
  }

  # Import all samples
  all_samples <- lapply(sample_files, import_oarfish_quant)

  # Combine into one data frame
  combined_counts <- do.call(cbind, all_samples)

  # Ensure column names are unique
  colnames(combined_counts) <- make.unique(colnames(combined_counts))

  # Add transcript_id as first column
  combined_counts$transcript_id <- rownames(combined_counts)
  combined_counts <- combined_counts[, c("transcript_id", setdiff(colnames(combined_counts), "transcript_id"))]

  # Sort by transcript_id
  combined_counts <- combined_counts %>% arrange(transcript_id)

  # Calculate proportions per sample
  tx_ids <- combined_counts$transcript_id
  counts_matrix <- combined_counts[, -1]
  proportions_df <- sweep(counts_matrix, 2, colSums(counts_matrix), "/")
  proportions_df$transcript_id <- tx_ids
  proportions_df <- proportions_df[, c("transcript_id", colnames(proportions_df)[-ncol(proportions_df)])]

  # Apply minimum proportion and sample thresholds
  filtered_transcripts <- proportions_df %>%
    filter(rowSums(select_if(., is.numeric) >= proportion_min) >= samples_min)

  # Subset counts accordingly
  combined_counts_filtered <- combined_counts %>% filter(transcript_id %in% filtered_transcripts$transcript_id)

  # Recalculate proportions on filtered data
  tx_ids <- combined_counts_filtered$transcript_id
  counts_matrix <- combined_counts_filtered[, -1]
  proportions_final <- sweep(counts_matrix, 2, colSums(counts_matrix), "/")
  proportions_final$transcript_id <- tx_ids
  proportions_final <- proportions_final[, c("transcript_id", colnames(proportions_final)[-ncol(proportions_final)])]

  # Calculate CPM
  CPM_df <- proportions_final
  CPM_df[, -1] <- CPM_df[, -1] * 1e6

  # Add gene_id column
  combined_counts_filtered <- cbind(gene_id = gene_id, combined_counts_filtered)
  proportions_final <- cbind(gene_id = gene_id, proportions_final)
  CPM_df <- cbind(gene_id = gene_id, CPM_df)

  # Export results
  write.csv(combined_counts_filtered, file.path(outdir, paste0(outdir, "_counts.csv")), quote = FALSE, row.names = FALSE)
  write.csv(proportions_final, file.path(outdir, paste0(outdir, "_proportions.csv")), quote = FALSE, row.names = FALSE)
  write.csv(CPM_df, file.path(outdir, paste0(outdir, "_CPM_values.csv")), quote = FALSE, row.names = FALSE)

  # Save total reads remaining
  total_reads_remaining <- sum(rowSums(select_if(combined_counts_filtered, is.numeric)))
  write.table(total_reads_remaining, file.path(outdir, "temp_files", "remaining_read_sum.txt"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)

})
