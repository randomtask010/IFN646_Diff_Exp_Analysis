library(ggplot2)
library(readxl)

file_path <- "Working Directory/Output/Summary_DE_Genes_by_Samples_Stats.xlsx"

samples <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")

# Iterate over each sample
for (sample in samples) {
  
  # Read the sheets for upregulated and downregulated genes
  upregulated <- read_excel(file_path, sheet = paste(sample, "fp_upregulated", sep = "_"))
  upregulated$regulation <- "Upregulated"
  
  downregulated <- read_excel(file_path, sheet = paste(sample, "fp_downregulated", sep = "_"))
  downregulated$regulation <- "Downregulated"
  
  # Combine the upregulated and downregulated data
  combined_data <- rbind(upregulated, downregulated)
  
  # Plotting
  p <- ggplot(combined_data, aes(x = gene_id)) +
    geom_point(aes(y = mean_condition_1, color = "Condition 1")) +  
    geom_point(aes(y = mean_condition_2, color = "Condition 2")) +  
    geom_line(aes(y = mean_condition_1, color = "Condition 1"), group = 1) +  
    geom_line(aes(y = mean_condition_2, color = "Condition 2"), group = 1) +  
    facet_wrap(~ regulation, scales = "free_y") +
    labs(
      title = paste("Mean Values in Sample", sample),
      x = "Gene ID",
      y = "Mean Expression Level"
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_color_manual(values = c("Condition 1" = "blue", "Condition 2" = "red"))
  
  
  
  print(p)
  
  output_file <- paste0("Working Directory/Output/Images/Mean_Diff_DownvsUp_FP_", sample, ".png")
  
 
  ggsave(output_file, plot = p, width = 10, height = 6)
}
