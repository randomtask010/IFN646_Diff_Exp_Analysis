run_loop_UpRegulation <- function(PValue, QValue) {
  
  library(VennDiagram)
  
  # Parameterised
  dir_path <- "Working Directory/Output/Threshold_Analysis/"
  tools <- c("deseq", "edgeR", "noiseq")
  samples <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")
  outlier_condition <- "outliers_upregulated"  
  output_image_dir <- "Working Directory/Output/Images/Threshold/"
  
  read_and_process_data <- function(tool, sample, PValue, QValue) {
    if (tool == "noiseq") {
      value_used <- QValue
      value_label <- "QValue"
    } else {
      value_used <- PValue
      value_label <- "PValue"
    }
    file_name <- paste0(dir_path, tool, "_", sample, "_", outlier_condition, "_", value_label, "_", value_used, ".csv")
    df <- read.csv(file_name, stringsAsFactors = FALSE)
    values <- df[[1]]
    return(values)
  }
  
  generate_and_save_venn <- function(sample, tool_values_list, PValue, QValue) {
    max_rows <- max(sapply(tool_values_list, length))
    
    # Extend each list to max_rows
    tool_values_list <- lapply(tool_values_list, function(values) {
      c(values, rep(NA, max_rows - length(values)))
    })
    
    filename = paste0(output_image_dir, "venn_", sample, "_", outlier_condition, "_PValue_", PValue, "_QValue_", QValue, ".png")
    
    # Create the Venn diagram
    venn.diagram(
      x = lapply(tool_values_list, function(values) values[!is.na(values)]),
      category.names = tools,
      output = TRUE,
      filename = filename,
      output.type = "png",
      imagetype = "png",
      resolution = 300,
      category.col = c("red", "blue", "green"),
      fill = c("red", "blue", "green")
    )
  }
  
  # Loop through samples
  for (sample in samples) {
    tool_values_list <- lapply(tools, function(tool) read_and_process_data(tool, sample, PValue, QValue))
    generate_and_save_venn(sample, tool_values_list, PValue, QValue)
  }
}
