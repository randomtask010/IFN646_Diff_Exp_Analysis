library(writexl)

# Load the original count data matrix for a sample
load_count_data <- function(sample) {
  count_data <- read.table(paste0("RAW Data/", sample, ".tsv"), header=TRUE, row.names=1)
  return(count_data)
}

# Load truth data for a sample
load_truth_data <- function(sample) {
  truth_data <- read.table(paste0("RAW Data/", sample, "_meta.tsv"), header=TRUE, row.names=1)
  return(truth_data)
}

samples <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")
conditions <- c("outliers_upregulated", "outliers_downregulated")

# List for summary stats
summary_stats <- list()

# Loop over each sample and condition
for (sample in samples) {
  
  # Load the truth data
  truth_data <- load_truth_data(sample)
  
  for (condition in conditions) {
    
    # Get relevant genes based on condition
    if (condition == "outliers_upregulated") {
      genes_list <- rownames(truth_data[truth_data$upregulation == 1,])
    } else if (condition == "outliers_downregulated") {
      genes_list <- rownames(truth_data[truth_data$downregulation == 1,])
    } else {
      next
    }
    
    # Load the original count data for the sample
    count_data <- load_count_data(sample)
    
    # Extract rows corresponding to the genes
    extracted_data <- count_data[rownames(count_data) %in% genes_list, ]
    
    # Determine sample size
    samplesize <- ncol(extracted_data) / 2
    
    # Calculate summary statistics for each gene
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
    
    # Store the statistics list with correct naming convention
    df <- as.data.frame(t(gene_stats))
    df$gene_id <- rownames(df)
    adjusted_condition <- ifelse(condition == "outliers_upregulated", "Truth_upregulated", 
                                 ifelse(condition == "outliers_downregulated", "Truth_downregulated", condition))
    summary_stats[[paste0(sample, "_", adjusted_condition)]] <- df
  }
}

# View the summary statistics
summary_stats
write_xlsx(summary_stats, "Working Directory/Output/Summary_DE_Genes_Truth_by_Samples_Stats.xlsx")
