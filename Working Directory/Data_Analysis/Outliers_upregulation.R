run_loop_UpRegulation <-function(PValue, QValue) {
  
  # Parameterized
  dir_path <- "Working Directory/Output/"
  tools <- c("deseq2", "edgeR", "noiseq")
  samples <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")
  outlier_condition <- "outliers_upregulated"  
  output_image_dir <- "Working Directory/Output/Images/"
  
  read_and_process_data <- function(tool, sample, PValue, QValue) {
    if (tool == "noiseq") {
      value_used <- QValue
    } else {
      value_used <- PValue
    }
    file_name <- paste0(dir_path, tool, "_", sample, "_", outlier_condition, "_PValue_", value_used, ".csv")
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
    
    if (tool == "noiseq") {
      value_used <- QValue
    } else {
      value_used <- PValue
    }
    filename = paste0(output_image_dir, "venn_", sample, "_", outlier_condition, "_PValue_", value_used, ".png")
    
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