library(ggplot2)

# Filter rows with Experiment values "3_500_500", "3_750_250", and "3_1000_0"
samplesize_3_df <- subset(Metrics_Final_Combined_DF, Experiment %in% c("3_500_500", "3_750_250", "3_1000_0"))

#01.Accuracy on sample size 3
p_acc_3 <- ggplot(samplesize_3_df, aes(x = Experiment, y = Accuracy, group = Tool, color = Tool)) +
  geom_line() +
  geom_point() +
  labs(title = "Gene count Impact on Accuracy Across Tools", x = "Gene count", y = "Accuracy") +
  theme_bw() +
  theme(legend.position = "top")  # Place the legend at the top
p_acc_3 <- p_acc_3 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_acc_3)

p_prec_3 <- ggplot(samplesize_3_df, aes(x = Experiment, y = Precision, group = Tool, color = Tool)) +
  geom_line() +
  geom_point() +
  labs(title = "Gene count Impact on Precision Across Tools", x = "Gene count", y = "Precision") +
  theme_bw() +
  theme(legend.position = "top")  # Place the legend at the top
p_prec_3 <- p_prec_3 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_prec_3)

#03.Recall == Sensitivity
p_rec_3 <- ggplot(samplesize_3_df, aes(x = Experiment, y = Recall, group = Tool, color = Tool)) +
  geom_line() +
  geom_point() +
  labs(title = "Gene count Impact on Recall(Sensitivity) Across Tools", x = "Gene count", y = "Recall") +
  theme_bw() +
  theme(legend.position = "top")  # Place the legend at the top
p_rec_3 <- p_rec_3 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_rec_3)

#04.F1_Score
p_f1_3 <- ggplot(samplesize_3_df, aes(x = Experiment, y = F1_Score, group = Tool, color = Tool)) +
  geom_line() +
  geom_point() +
  labs(title = "Gene count Impact on F1_Score Across Tools", x = "Gene count", y = "F1_Score") +
  theme_bw() +
  theme(legend.position = "top")  # Place the legend at the top
p_f1_3 <- p_f1_3 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p_f1_3)

ggsave(filename = "Working Directory/Output/Images/sample3_GeneCountVsAccuracy_Plot.png", plot = p_acc_3, width = 6, height = 6)
ggsave(filename = "Working Directory/Output/Images/sample3_GeneCountVsPrecision_Plot.png", plot = p_prec_3, width = 6, height = 6)
ggsave(filename = "Working Directory/Output/Images/sample3_GeneCountVsRecall_Plot.png", plot = p_rec_3, width = 6, height = 6)
ggsave(filename = "Working Directory/Output/Images/sample3_GeneCountVsF1_Score_Plot.png", plot = p_f1_3, width = 6, height = 6)

