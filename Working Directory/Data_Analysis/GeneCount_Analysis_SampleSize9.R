library(ggplot2)

# Filter rows with Experiment values "9_500_500", "9_750_250", and "9_1000_0"
samplesize_9_df <- subset(Metrics_Final_Combined_DF, Experiment %in% c("9_500_500", "9_750_250", "9_1000_0"))

#01.Accuracy on sample size 9
p_acc_9 <- ggplot(samplesize_9_df, aes(x = Experiment, y = Accuracy, group = Tool, color = Tool)) +
  geom_line() +
  geom_point() +
  labs(title = "Gene count Impact on Accuracy Across Tools", x = "Gene count", y = "Accuracy") +
  theme_bw() +
  theme(legend.position = "top")  # Place the legend at the top
p_acc_9 <- p_acc_9 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_acc_9)

#Precision
p_prec_9 <- ggplot(samplesize_9_df, aes(x = Experiment, y = Precision, group = Tool, color = Tool)) +
  geom_line() +
  geom_point() +
  labs(title = "Gene count Impact on Precision Across Tools", x = "Gene count", y = "Precision") +
  theme_bw() +
  theme(legend.position = "top")  # Place the legend at the top
p_prec_9 <- p_prec_9 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_prec_9)

#03.Recall == Sensitivity
p_rec_9 <- ggplot(samplesize_9_df, aes(x = Experiment, y = Recall, group = Tool, color = Tool)) +
  geom_line() +
  geom_point() +
  labs(title = "Gene count Impact on Recall(Sensitivity) Across Tools", x = "Gene count", y = "Recall") +
  theme_bw() +
  theme(legend.position = "top")  # Place the legend at the top
p_rec_9 <- p_rec_9 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_rec_9)

#04.F1_Score
p_f1_9 <- ggplot(samplesize_9_df, aes(x = Experiment, y = F1_Score, group = Tool, color = Tool)) +
  geom_line() +
  geom_point() +
  labs(title = "Gene count Impact on F1_Score Across Tools", x = "Gene count", y = "F1_Score") +
  theme_bw() +
  theme(legend.position = "top")  # Place the legend at the top
p_f1_9 <- p_f1_9 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_f1_9)

ggsave(filename = "Working Directory/Output/Images/sample9_GeneCountVsAccuracy_Plot.png", plot = p_acc_9, width = 6, height = 6)
ggsave(filename = "Working Directory/Output/Images/sample9_GeneCountVsPrecision_Plot.png", plot = p_prec_9, width = 6, height = 6)
ggsave(filename = "Working Directory/Output/Images/sample9_GeneCountVsRecall_Plot.png", plot = p_rec_9, width = 6, height = 6)
ggsave(filename = "Working Directory/Output/Images/sample9_GeneCountVsF1_Score_Plot.png", plot = p_f1_9, width = 6, height = 6)

