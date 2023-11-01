library(ggplot2)
library(readxl)

file_path <- "Working Directory/Output/Summary_DE_Genes_Truth_by_Samples_Stats.xlsx"
samples <- c("3_500_500", "3_750_250", "6_500_500", "6_750_250", "9_500_500", "9_750_250")

# Iterate over each sample
for (sample in samples) {
  
  # Read the sheets for upregulated and downregulated genes
  upregulated <- read_excel(file_path, sheet = paste(sample, "Truth_upregulated", sep = "_"))
  upregulated$regulation <- "Upregulated"
  
  downregulated <- read_excel(file_path, sheet = paste(sample, "Truth_downregulated", sep = "_"))
  downregulated$regulation <- "Downregulated"
  
  # Combine the upregulated and downregulated data
  combined_data <- rbind(upregulated, downregulated)
  
  # Calculate M and A values
  combined_data$M <- log2(combined_data$mean_condition_1 / combined_data$mean_condition_2)
  combined_data$A <- 0.5 * (log2(combined_data$mean_condition_1) + log2(combined_data$mean_condition_2))
  
  # Plotting MA-plot
  p <- ggplot(combined_data, aes(x = A, y = M)) +
    geom_point(aes(color = regulation), alpha = 0.7) +
    facet_wrap(~ regulation, scales = "free_y") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
    labs(
      title = paste("MA-Plot for Sample", sample),
      x = "A (Average Expression)",
      y = "M (Log2 Fold Change)"
    ) +
    theme_bw() +
    scale_color_manual(values = c("Upregulated" = "blue", "Downregulated" = "red"))
  
  print(p)
  
  output_file <- paste0("Working Directory/Output/Images/MA_Plot_Truth", sample, ".png")
  ggsave(output_file, plot = p, width = 10, height = 6)
}
