## Incorporate sample size into this!

# Function to load the original count data matrix for a sample
load_count_data <- function(sample) {
  count_data <- read.table(paste0("RAW Data/", sample, ".tsv"), header=TRUE, row.names=1)
  return(count_data)
}

samples <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")
conditions <- c("outliers_upregulated", "outliers_downregulated")

# List for  summary stats
summary_stats <- list()

# Loop over each sample and condition
for (sample in samples) {
  for (condition in conditions) {
    
    # Load the genes from the generated CSV
    file_path <- paste0("Working Directory/Output/Common_DE_Genes_Thresholds_", condition, "_", sample, ".csv")
    if (!file.exists(file_path)) next
    genes_df <- read.csv(file_path, stringsAsFactors = FALSE)
    
    # Convert the dataframe to a list of genes
    genes_list <- unlist(genes_df, use.names = FALSE)
    genes_list <- genes_list[!is.na(genes_list)]
    
    # Load the original count data for the sample
    count_data <- load_count_data(sample)
    
    # Extract rows corresponding to the genes
    extracted_data <- count_data[rownames(count_data) %in% genes_list, ]
    
    # Calculate summary statistics for each gene
    gene_stats <- apply(extracted_data, 1, function(gene_counts) {
      c(
        mean = mean(gene_counts, na.rm=TRUE),
        median = median(gene_counts, na.rm=TRUE),
        sd = sd(gene_counts, na.rm=TRUE),
        min = min(gene_counts, na.rm=TRUE),
        Q1 = quantile(gene_counts, 0.25, na.rm=TRUE),
        Q3 = quantile(gene_counts, 0.75, na.rm=TRUE),
        max = max(gene_counts, na.rm=TRUE)
      )
    })
    
    # Store the statistics in the list
    summary_stats[[paste0(sample, "_", condition)]] <- t(gene_stats)
  }
}

# View the summary statistics
summary_stats
