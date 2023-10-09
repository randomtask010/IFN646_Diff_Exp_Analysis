library(ggplot2)

ggplot(Metrics_Final_Combined_DF, aes(x = Experiment, y = Accuracy, fill = Tool)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Accuracy of Tools Across Experiments")
