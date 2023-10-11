
library(ggplot2)

# plot setup
p <- ggplot(Metrics_Final_Combined_DF, aes(x = Experiment, y = Accuracy, fill = Tool))
p <- p + geom_bar(stat = "identity", position = "dodge")
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 45, hjust = 1))
p <- p + labs(title = "Accuracy of Tools Across Experiments")


# Print the plot 
print(p)
#Save to Output/Images folder
ggsave(filename = "Working Directory/Output/Images/Accuracy_Metrics_Plot.png", plot = p)
