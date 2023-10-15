library(ggplot2)
library(dplyr)

# Load the data
edger_data <- read.csv("Working Directory/Output/Threshold_Analysis/Threshold_edger.csv")
noiseq_data <- read.csv("Working Directory/Output/Threshold_Analysis/Threshold_noiseq.csv")
deseq_data <- read.csv("Working Directory/Output/Threshold_Analysis/Threshold_deseq.csv")

# Add a 'Tool' column to each dataset
edger_data$Tool <- 'edgeR'
noiseq_data$Tool <- 'NOISeq'
deseq_data$Tool <- 'DESeq2'

# Combine the datasets
combined_data <- bind_rows(edger_data, noiseq_data, deseq_data)

# Plot
ggplot(combined_data, aes(x=Threshold, y=FDR, color=Tool)) +
  geom_line() +
  facet_wrap(~Experiment, scales="free_x") +
  labs(title="Comparison of FDR across Tools at different Thresholds",
       x="Threshold",
       y="False Discovery Rate (FDR)") +
  theme_minimal()
