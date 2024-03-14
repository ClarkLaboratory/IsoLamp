#!/usr/bin/env Rscript


suppressPackageStartupMessages({
  library(FactoMineR)
  library(factoextra)
  library(ggplot2)
  library(ggrepel)
 library(optparse) 
})


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
  
  count_data <- (opt$input)
  output_prefix <- (opt$output_prefix)
  grouping_variable <- (opt$grouping_variable)


  # Check if the grouping_variable option was provided and is not "NULL"
  if (!is.null(opt$grouping_variable) && opt$grouping_variable != "NULL") {
    grouping_variable <- opt$grouping_variable
  } else {
    grouping_variable <- NULL
  }

plot_pca <- function(count_data, output_prefix, grouping_variable) {
  
  # Read in the counts table from a CSV file
  data <- read.csv(count_data, header = TRUE, row.names = 2)
  
  #transpose the data
  data_transposed <- t(data[,-1])
  
  # Perform PCA on the transposed data
  pca_result <- PCA(data_transposed, graph = FALSE)
  
  # Create a data frame for PCA results
  pca_data <- as.data.frame(pca_result$ind$coord)
  
  # Rename columns for the PCA data frame (PC1, PC2, etc.)
  colnames(pca_data) <- c(paste0("PC", 1:ncol(pca_data)))
  pca_data$sample <- rownames(pca_data)
  
  if (is.null(grouping_variable)) {
    # If grouping_variable is NULL, set the "Group" column to 1 for all samples
    pca_data$group <- "Group A"
    merged_df <- pca_data
  } else {
    # Read in the grouping variable
    d.f <- read.csv(grouping_variable, header = TRUE)
    merged_df <- merge(d.f, pca_data, by = "sample", all = TRUE)
  }
 
  # Create a PCA plot using ggplot2
  pca_plot <- ggplot(merged_df, aes(x = PC1, y = PC2, color = group, label = merged_df$sample)) +
    geom_point(size = 3, alpha = 0.6) +  # Customize point size and color
    geom_text_repel(
      hjust = -0.2, vjust = -0.5, size = 4, color = "black",
      segment.color = "transparent"  # Make the connecting lines transparent
    ) +
    labs(
      caption = paste("Source:", count_data) # Add a caption
    ) +
    xlab(paste0("PC1: ", round(pca_result[["eig"]][1, 2]), "% variance")) +
    ylab(paste0("PC2: ", round(pca_result[["eig"]][2, 2]), "% variance")) +
    theme_bw() +  # Use a minimal theme
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),  # Adjust title appearance
      plot.caption = element_text(hjust = 1, size = 10)  # Adjust caption appearance
    )

  # Return the PCA plot
  pdf(paste0(output_prefix, "_PCA.pdf"), width = 8, height = 8)
  plot(pca_plot)
  dev.off()

}

plot_pca(count_data, output_prefix, grouping_variable)