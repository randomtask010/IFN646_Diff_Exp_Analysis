# List of file paths
file_paths <- c(
  "Working Directory/Output/Threshold_noiseq_3_1000_0.csv",
  "Working Directory/Output/Threshold_noiseq_6_1000_0.csv",
  "Working Directory/Output/Threshold_noiseq_9_1000_0.csv"
)

# Initialize an empty dataframe to store the merged data
merged_data <- data.frame()

# Loop through each file and merge its data into the merged_data dataframe
for (file_path in file_paths) {
  # Read the CSV file into a temporary dataframe
  temp_data <- read.csv(file_path, header = TRUE)
  
  # Combine the temporary dataframe with the merged_data dataframe
  merged_data <- rbind(merged_data, temp_data)
}
merged_data
# Now, merged_data contains the combined data from all three CSV files


p <- ggplot(merged_data, aes(x = Threshold, y = FDR, group = Experiment, color = Experiment)) +
  geom_line() +
  geom_point() +
  labs(title = "q-value threshold against FDR with NOISeq tool", x = "q-value threshold", y = "FDR") +
  theme_bw() 
#theme(legend.position = "top")  # Place the legend at the top

# Rotate x-axis labels
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
print(p)
ggsave(filename = "Working Directory/Output/Images/Threshold_FDR_Noiseq_Plot_3_6_9_1000_0.png", plot = p, width = 6, height = 6)



p_rec <- ggplot(merged_data, aes(x = Threshold, y = Recall, group = Experiment, color = Experiment)) +
  geom_line() +
  geom_point() +
  labs(title = "q-value threshold against Recall with NOISeq tool", x = "q-value threshold", y = "Recall") +
  theme_bw() 
#theme(legend.position = "top")  # Place the legend at the top

# Rotate x-axis labels
p_rec <- p_rec + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
print(p_rec)
ggsave(filename = "Working Directory/Output/Images/Threshold_Recall_Noiseq_Plot_3_6_9_1000_0.png", plot = p_rec, width = 6, height = 6)
