library(writexl)

# Load the original count data matrix for a sample
load_count_data <- function(sample) {
  count_data <- read.table(paste0("RAW Data/", sample, ".tsv"), header=TRUE, row.names=1)
  return(count_data)
}

samples <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")
conditions <- c("outliers_upregulated", "outliers_downregulated")

# List for summary stats
summary_stats <- list()

# Loop over each sample and condition
for (sample in samples) {
  for (condition in conditions) {
    
    # Load the genes from the CSVs
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
    
    # Determine sample size
    samplesize <- ncol(extracted_data) / 2
    
    # Calculate summary statistics for each gene divided by sample size
    # This allows us to identify if the samples potentially were so skewed that normalisation didnt account for them leading to consistently DE genes
    gene_stats <- apply(extracted_data, 1, function(gene_counts) {
      current_sample <- gene_counts[1:samplesize]
      next_sample <- gene_counts[(samplesize+1):(2*samplesize)]
      c(
        mean_condition_1 = mean(current_sample, na.rm=TRUE),
        median_condition_1 = median(current_sample, na.rm=TRUE),
        sd_condition_1 = sd(current_sample, na.rm=TRUE),
        min_condition_1 = min(current_sample, na.rm=TRUE),
        Q1_condition_1 = quantile(current_sample, 0.25, na.rm=TRUE),
        Q3_condition_1 = quantile(current_sample, 0.75, na.rm=TRUE),
        max_condition_1 = max(current_sample, na.rm=TRUE),
        mean_condition_2 = mean(next_sample, na.rm=TRUE),
        median_condition_2 = median(next_sample, na.rm=TRUE),
        sd_condition_2 = sd(next_sample, na.rm=TRUE),
        min_condition_2 = min(next_sample, na.rm=TRUE),
        Q1_condition_2 = quantile(next_sample, 0.25, na.rm=TRUE),
        Q3_condition_2 = quantile(next_sample, 0.75, na.rm=TRUE),
        max_condition_2 = max(next_sample, na.rm=TRUE)
      )
    })
    
    # Store the statistics list
    df <- as.data.frame(t(gene_stats))
    df$gene_id <- rownames(df)
    summary_stats[[paste0(sample, "_", condition)]] <- df
    
  }
}

# View the summary statistics
summary_stats
write_xlsx(summary_stats, "Working Directory/Output/Summary_DE_Genes_by_Samples_Stats.xlsx")