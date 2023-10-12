library(ggplot2)

# plot setup
p <- ggplot(Metrics_Final_Combined_DF, aes(x = Experiment, y = Accuracy, fill = Tool))
p <- p + geom_bar(stat = "identity", position = "dodge")
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p <- p + labs(title = "Accuracy of Tools Across Experiments", x = "Experiment", y = "Accuracy")

# plot setup
p_fdr <- ggplot(Metrics_Final_Combined_DF, aes(x = Experiment, y = FDR, fill = Tool))
p_fdr <- p_fdr + geom_bar(stat = "identity", position = "dodge")
p_fdr <- p_fdr + theme_bw()
p_fdr <- p_fdr + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_fdr <- p_fdr + labs(title = "False Discovery Rate of Tools Across Experiments",x = "Experiment", y = "False Discovery Rate")

#Sensitivity
p_rec <- ggplot(Metrics_Final_Combined_DF, aes(x = Experiment, y = Recall, fill = Tool))
p_rec <- p_rec + geom_bar(stat = "identity", position = "dodge")
p_rec <- p_rec + theme_bw()
p_rec <- p_rec + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_rec <- p_rec + labs(title = "Sensitivity of Tools Across Experiments", x = "Experiment", y = "Sensitivity")



#Save to Output/Images folder
ggsave(filename = "Working Directory/Output/Images/Accuracy_Metrics_Plot.png", plot = p, width = 6, height = 6)
ggsave(filename = "Working Directory/Output/Images/FDR_Metrics_Plot.png", plot = p_fdr, width = 6, height = 6)
ggsave(filename = "Working Directory/Output/Images/Recall_Metrics_Plot.png", plot = p_rec, width = 6, height = 6)

