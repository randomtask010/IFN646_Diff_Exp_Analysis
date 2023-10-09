
# paramertised
dir_path <- "Working Directory/Output/"
tools <- c("deseq2", "edgeR", "noiseq")
samples <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")

data <- data.frame()

# Loop for parameters
for (tool in tools) {
  for (sample in samples) {
    file_name <- paste0(dir_path, "Metrics_", tool, "_", sample, ".csv")
    # Read the CSV
    df <- read.csv(file_name, stringsAsFactors = FALSE)
    df$Tool <- tool
    df$Experiment <- sample
    # Stack the data verticlaly
    data <- rbind(data, df)
  }
}

#Rename for useage and analysis
Metrics_Final_Combined_DF <- data



