# Load the required package
library(ggplot2)

# Replace 'your_file.xlsx' with the actual file path
file_path <- "Working Directory/Output/Summary_DE_Genes_by_Samples_Stats.xlsx"

# Get the sheet names
sheet_names <- excel_sheets(file_path)  # For 'readxl' package

# Create a list to store the data frames
data_list <- list()

## Draw mean value plots 

# Read each sheet into a data frame and store it in the list
for (sheet in sheet_names) {
  data_list[[sheet]] <- read_excel(file_path, sheet = sheet)  
  
  # Create a bar plot to visualize the difference between means for each sheet
  p <- ggplot(data_list[[sheet]], aes(x = gene_id, y = mean_condition_2 - mean_condition_1)) +
    geom_bar(stat = "identity", fill = "blue") +
    labs(
      title = paste("Difference between Mean Condition 2 and Mean Condition 1 for", sheet),
      x = "Gene ID",
      y = "Mean Condition 2 - Mean Condition 1"
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(labels = scales::comma) +
    coord_cartesian(ylim = c(-100, 100))  # Adjust the y-axis limits to zoom in on small differences
  
  # Print or save the plot for each sheet
  output_directory <- "Working Directory/Output/Images/"
  output_file <- file.path(output_directory, paste("MeanDiff_", sheet, ".png"))
  ggsave(output_file, plot = p, width = 10, height = 4)  
  
  # Print the saved file path
  print(output_file)
}



###Draw SD value plots

# Read each sheet into a data frame and store it in the list
for (sheet in sheet_names) {
  data_list[[sheet]] <- read_excel(file_path, sheet = sheet)
  
  # Create a bar plot to visualize the difference between means for each sheet
  p <- ggplot(data_list[[sheet]], aes(x = gene_id, y = sd_condition_2 - sd_condition_1)) +
    geom_bar(stat = "identity", fill = "blue") +
    labs(
      title = paste("Difference between sd Condition 2 and sd Condition 1 for", sheet),
      x = "Gene ID",
      y = "sd Condition 2 - sd Condition 1"
    ) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    scale_y_continuous(labels = scales::comma) +
    coord_cartesian(ylim = c(-100, 100))  # Adjust the y-axis limits to zoom in on small differences
  
  output_directory <- "Working Directory/Output/Images/"
  output_file <- file.path(output_directory, paste("SdDiff_", sheet, ".png"))
  ggsave(output_file, plot = p, width = 10, height = 4) 
  
  # Print the saved file path
  print(output_file)
}