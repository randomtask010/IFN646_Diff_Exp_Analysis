library(ggplot2)

# List of file paths

threshold_edger_data <- read.csv("Working Directory/Output/Threshold_edger.csv", header = TRUE)

#q-value threshold against FDR
p_fdr <- ggplot(threshold_edger_data, aes(x = Threshold, y = FDR, group = Experiment, color = Experiment)) +
  geom_line() +
  geom_point() +
  labs(title = "p-value threshold against FDR with Edger tool", x = "p-value threshold", y = "False Discovery Rate") +
  theme_bw() 
ggsave(filename = "Working Directory/Output/Images/Threshold_FDR_Edger_Plot.png", plot = p_fdr, width = 6, height = 6)

#q-value threshold against Recall(Sensitivity)
p_rec <- ggplot(threshold_edger_data, aes(x = Threshold, y = Recall, group = Experiment, color = Experiment)) +
  geom_line() +
  geom_point() +
  labs(title = "p-value threshold against Recall(Sensitivity) with Edger tool", x = "p-value threshold", y = "Sensitivity") +
  theme_bw() 
ggsave(filename = "Working Directory/Output/Images/Threshold_Recall_Edger_Plot.png", plot = p_rec, width = 6, height = 6)

#q-value threshold against Accuracy
p_acc <- ggplot(threshold_edger_data, aes(x = Threshold, y = Accuracy, group = Experiment, color = Experiment)) +
  geom_line() +
  geom_point() +
  labs(title = "p-value threshold against Accuracy with Edger tool", x = "p-value threshold", y = "Accuracy") +
  theme_bw() 
ggsave(filename = "Working Directory/Output/Images/Threshold_Accuracy_Edger_Plot.png", plot = p_acc, width = 6, height = 6)


#Check sensitivity of sample size vs threshold impact
filtered_df_500_500 <- subset(threshold_edger_data, Experiment %in% c("3_500_500", "6_500_500", "9_500_500"))
p_sen <- ggplot(filtered_df_500_500, aes(x = Threshold, y = Recall, group = Experiment, color = Experiment)) +
  geom_line() +
  geom_point() +
  labs(title = "p-value threshold against Recall(Sensitivity) with Edger tool", x = "p-value threshold", y = "Sensitivity") +
  theme_bw() 
ggsave(filename = "Working Directory/Output/Images/Threshold_Recall_Edger_Plot_3_6_9_500_500.png", plot = p_sen, width = 6, height = 6)


#q-value threshold against Recall bar graph
p_bar <- ggplot(threshold_edger_data, aes(x = Threshold, y = Recall, fill = Experiment)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_bw() +
  labs(title = "p-value threshold against Recall(Sensitivity) with Edger tool", x = "p-value threshold", y = "Sensitivity")
ggsave(filename = "Working Directory/Output/Images/Threshold_Recall_Edger_Bar_Plot.png", plot = p_bar, width = 10, height = 10)
