#!/usr/bin/env Rscript

# in IsoLamp main, if grouping_variable is NULL, don't run Rscript

suppressPackageStartupMessages({
  library(reshape)
  library(dplyr)
  library(rstatix)
  library(purrr)
  library(optparse) 
})

# importing
option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="input count_data"),
  make_option(c("-o", "--output_prefix"), type="character", default=NULL,
              help="output name"),
  make_option(c("-g", "--grouping_variable"), type = "character", default = NULL,
              help = "grouping table")
)

opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

count_props <- (opt$input)
output_prefix <- (opt$output_prefix)
grouping_variable <- (opt$grouping_variable)

groups <- suppressWarnings(read.csv(grouping_variable, header=T))
counts <- suppressWarnings(read.csv(count_props, header=T))

# melt df into long format
counts_melt <- reshape::melt(counts, variable_name="sample")
# merge with groups variable
counts_merged <- merge(counts_melt, groups, by="sample", all.x=T)

# function to perform t test
compare_proportions_between_groups <- function(df) {
  #df <- list_tx[[1]]
  # store TXNAME
  tx_id <- df$transcript_id
  
  # run t test
  tmp <- tryCatch(
        {
          df %>% rstatix::t_test(value ~ group, p.adjust.method = "fdr", detailed = TRUE, paired = FALSE, alternative = "two.sided")
        },
        error = function(e) {
          return(NA)
        }
      )
  
  # catch NA from tmp if error
  if (all(is.na(tmp))) {
        return(NA)
  } else {
        tmp
        result <- tmp %>% 
          dplyr::mutate(transcript_id = tx_id[1]) %>% 
          subset(select = c("transcript_id", "group1", "group2", "estimate", "estimate1", "estimate2", "n1", "n2", "p")) %>% 
          dplyr::rename("mean_prop_diff" = "estimate",
                        "group1_mean_prop" = "estimate1",
                        "group2_mean_prop" = "estimate2",
                        "group1_obs" = "n1",
                        "group2_obs" = "n2",
                        "padj" = "p")
  }
  
  return(result)
}

# split counts df by TXNAME
list_tx <- split(counts_merged, counts_merged$transcript_id)

# apply t test function
results_list_tx <- purrr::map(list_tx, compare_proportions_between_groups)

if (all(is.na(results_list_tx))) {
  # if all are NA, don't output result
  print("No isoform proportions could be compared")
} else {
  # remove NA elements
  results_list_tx <- results_list_tx[!is.na(results_list_tx)]
  
  # combine results output into a dataframe
  results_list_tx_df <- bind_rows(results_list_tx, .id = c("transcript_id"))
  
  # output
  write.csv(results_list_tx_df, paste0(output_prefix, "/", output_prefix, "_prop_t_test_results.csv"), row.names=F, quote=F)
  print("Complete")
}


