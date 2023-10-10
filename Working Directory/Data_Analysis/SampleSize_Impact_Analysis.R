library(ggplot2)

# Filter rows with Experiment values "3_500_500", "6_500_500", and "9_500_500"
filtered_df <- subset(Metrics_Final_Combined_DF, Experiment %in% c("3_500_500", "6_500_500", "9_500_500"))
filtered_df

#01.Accuracy
p_acc <- ggplot(filtered_df, aes(x = Experiment, y = Accuracy, group = Tool, color = Tool)) +
  geom_line() +
  geom_point() +
  labs(title = "Sample Size Impact on Accuracy Across Tools", x = "Sample Size", y = "Accuracy") +
  theme_bw() +
  theme(legend.position = "top")  # Place the legend at the top

# Rotate x-axis labels
p_acc <- p_acc + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
print(p_acc)


#02.Precision
p_prec <- ggplot(filtered_df, aes(x = Experiment, y = Precision, group = Tool, color = Tool)) +
  geom_line() +
  geom_point() +
  labs(title = "Sample Size Impact on Precision Across Tools", x = "Sample Size", y = "Precision") +
  theme_bw() +
  theme(legend.position = "top")  # Place the legend at the top

# Rotate x-axis labels
p_prec <- p_prec + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
print(p_prec)

#03.Recall == Sensitivity
p_rec <- ggplot(filtered_df, aes(x = Experiment, y = Recall, group = Tool, color = Tool)) +
  geom_line() +
  geom_point() +
  labs(title = "Sample Size Impact on Recall(Sensitivity) Across Tools", x = "Sample Size", y = "Recall") +
  theme_bw() +
  theme(legend.position = "top")  # Place the legend at the top

# Rotate x-axis labels
p_rec <- p_rec + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
print(p_rec)

#04.F1_Score
p_f1 <- ggplot(filtered_df, aes(x = Experiment, y = F1_Score, group = Tool, color = Tool)) +
  geom_line() +
  geom_point() +
  labs(title = "Sample Size Impact on F1_Score Across Tools", x = "Sample Size", y = "F1_Score") +
  theme_bw() +
  theme(legend.position = "top")  # Place the legend at the top

# Rotate x-axis labels
p_f1 <- p_f1 + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Print the plot
print(p_f1)

#Save plots to Output/Images folder
ggsave(filename = "Working Directory/Output/Images/SampleSizeVsAccuracy_Plot.png", plot = p_acc, width = 6, height = 6)
ggsave(filename = "Working Directory/Output/Images/SampleSizeVsPrecision_Plot.png", plot = p_prec, width = 6, height = 6)
ggsave(filename = "Working Directory/Output/Images/SampleSizeVsRecall_Plot.png", plot = p_rec, width = 6, height = 6)
ggsave(filename = "Working Directory/Output/Images/SampleSizeVsF1_Score_Plot.png", plot = p_f1, width = 6, height = 6)
