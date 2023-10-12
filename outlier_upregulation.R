# Parameterized
dir_path <- "Working Directory/Output/"
tools <- c("deseq2", "edgeR", "noiseq")
samples <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")
outlier_condition <- "outliers_upregulated"  
output_image_dir <- "Working Directory/Output/Images/"

generate_and_save_venn <- function(tool, sample) {
  file_name <- paste0(dir_path, tool, "_", sample, "_", outlier_condition, ".csv")
  df <- read.csv(file_name, stringsAsFactors = FALSE)
  values <- df[[1]]
  max_rows <- max(length(values))
  values <- c(values, rep(NA, max_rows - length(values)))
  
  data <- data.frame(matrix(nrow = length(values), ncol = 0))
  col_name <- paste0(tool, "_", sample, "_UPREGULATED")
  data[col_name] <- values
  
  # Create the Venn diagram
  venn.diagram(
    x = list(tool = values[!is.na(values)]),
    category.names = tool,
    output = TRUE,
    filename = paste0(output_image_dir, "venn_", sample, "_", tool, "_", outlier_condition, ".png"),
    output.type = "png",
    imagetype = "png",
    resolution = 300,
    category.col = c("red", "blue", "green"),
  )
  
  return(data)
}

# Loop through tools and samples
for (tool in tools) {
  for (sample in samples) {
    generate_and_save_venn(tool, sample)
  }
}

