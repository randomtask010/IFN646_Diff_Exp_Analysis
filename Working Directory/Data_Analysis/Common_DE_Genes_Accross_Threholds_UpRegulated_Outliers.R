# Function to gather DE genes from a specific file
gather_DE_genes <- function(tool, sample, value, type) {
  if (tool == "noiseq") {
    value_label <- "QValue"
  } else {
    value_label <- "PValue"
  }
  file_name <- paste0(dir_path, tool, "_", sample, "_outliers_upregulated_", value_label, "_", value, ".csv")
  if (!file.exists(file_name)) return(character(0))
  df <- read.csv(file_name, stringsAsFactors = FALSE)
  return(df[[1]])
}

dir_path <- "Working Directory/Output/Threshold_Analysis/"
tools <- c("deseq", "edgeR", "noiseq")
samples <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")
PValues <- c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09)
QValues <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)

common_genes_upregulated_by_threshold <- list()

for (i in 1:length(PValues)) {
  p <- PValues[i]
  q <- QValues[i]
  
  all_genes_by_tool <- list()
  for (tool in tools) {
    genes <- unlist(lapply(samples, gather_DE_genes, tool=tool, value=ifelse(tool == "noiseq", q, p), type=ifelse(tool == "noiseq", "QValue", "PValue")))
    all_genes_by_tool[[tool]] <- unique(genes)
  }
  
  common_genes <- Reduce(intersect, all_genes_by_tool)
  common_genes_upregulated_by_threshold[[paste0("PValue_", p, "_QValue_", q)]] <- common_genes
}

common_genes_upregulated_by_threshold

# Convert the list to a data frame
max_genes <- max(sapply(common_genes_upregulated_by_threshold, length))
df_upregulated <- as.data.frame(do.call(cbind, lapply(common_genes_upregulated_by_threshold, function(genes) {
  c(genes, rep(NA, max_genes - length(genes)))
})))

# Save the data frame to CSV
write.csv(df_upregulated, "Working Directory/Output/Common_DE_Genes_Thresholds_UpRegulated_Outliers.csv", row.names = FALSE)

