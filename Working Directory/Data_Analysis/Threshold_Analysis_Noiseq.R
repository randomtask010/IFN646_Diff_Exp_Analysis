library(ggplot2)

# List of file paths

threshold_noiseq_data <- read.csv("Working Directory/Output/Threshold_noiseq.csv", header = TRUE)

#q-value threshold against FDR
p_fdr <- ggplot(threshold_noiseq_data, aes(x = Threshold, y = FDR, group = Experiment, color = Experiment)) +
  geom_line() +
  geom_point() +
  labs(title = "q-value threshold against FDR with NOISeq tool", x = "q-value", y = "False Discovery Rate") +
  theme_bw() 
ggsave(filename = "Working Directory/Output/Images/Threshold_FDR_Noiseq_Plot.png", plot = p_fdr, width = 6, height = 6)

#q-value threshold against Recall(Sensitivity)
p_rec <- ggplot(threshold_noiseq_data, aes(x = Threshold, y = Recall, group = Experiment, color = Experiment)) +
  geom_line() +
  geom_point() +
  labs(title = "q-value threshold against Recall(Sensitivity) with NOISeq tool", x = "q-value threshold", y = "Sensitivity") +
  theme_bw() 
ggsave(filename = "Working Directory/Output/Images/Threshold_Recall_Noiseq_Plot.png", plot = p_rec, width = 6, height = 6)

#q-value threshold against Accuracy
p_acc <- ggplot(threshold_noiseq_data, aes(x = Threshold, y = Accuracy, group = Experiment, color = Experiment)) +
  geom_line() +
  geom_point() +
  labs(title = "q-value threshold against Accuracy with NOISeq tool", x = "q-value threshold", y = "Accuracy") +
  theme_bw() 
ggsave(filename = "Working Directory/Output/Images/Threshold_Accuracy_Noiseq_Plot.png", plot = p_acc, width = 6, height = 6)


#Check sensitivity of sample size vs threshold impact
filtered_df_500_500 <- subset(threshold_noiseq_data, Experiment %in% c("3_500_500", "6_500_500", "9_500_500"))
p_sen <- ggplot(filtered_df_500_500, aes(x = Threshold, y = Recall, group = Experiment, color = Experiment)) +
  geom_line() +
  geom_point() +
  labs(title = "q-value threshold against Recall(Sensitivity) with NOISeq tool", x = "q-value threshold", y = "Sensitivity") +
  theme_bw() 
ggsave(filename = "Working Directory/Output/Images/Threshold_Recall_Noiseq_Plot_3_6_9_500_500.png", plot = p_sen, width = 6, height = 6)


#q-value threshold against Recall bar graph
p_bar <- ggplot(threshold_noiseq_data, aes(x = Threshold, y = Recall, fill = Experiment)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(title = "q-value threshold against Recall(Sensitivity) with NOISeq tool", x = "q-value threshold", y = "Sensitivity")
ggsave(filename = "Working Directory/Output/Images/Threshold_Recall_Noiseq_Bar_Plot.png", plot = p_bar, width = 10, height = 10)
