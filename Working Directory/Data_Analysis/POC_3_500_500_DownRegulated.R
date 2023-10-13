# POC for 3_500_500 - THIS WILL NOT STAY, but its a 1am POC that will need to be paramertised and made elegant - we will eventually move away from the manual csv comparisons for up/down regulated comparisons of overlaping mis-identified genes
desq2_3_500_500_DOWN_df <- read.csv("Working Directory/Output/deseq2_baseline_3_500_500_outliers_downregulated_PValue_0.05.csv", stringsAsFactors = FALSE)
edger_3_500_500_DOWN_df <- read.csv("Working Directory/Output/edger_baseline_3_500_500_outliers_downregulated_PValue_0.05.csv", stringsAsFactors = FALSE)
noiseq_3_500_500_DOWN_df <- read.csv("Working Directory/Output/noiseq_baseline_3_500_500_outliers_downregulated_PValue_0.8.csv", stringsAsFactors = FALSE)
deseq2_values <- desq2_3_500_500_DOWN_df[[1]]
edger_values <- edger_3_500_500_DOWN_df[[1]]
noiseq_values <- noiseq_3_500_500_DOWN_df[[1]]
max_rows <- max(length(deseq2_values), length(edger_values), length(noiseq_values))
deseq2_values <- c(deseq2_values, rep(NA, max_rows - length(deseq2_values)))
edger_values <- c(edger_values, rep(NA, max_rows - length(edger_values)))
noiseq_values <- c(noiseq_values, rep(NA, max_rows - length(noiseq_values)))
downregulated_genes_3_500_500 <- data.frame(desq2_3_500_500_DOWNREGULATED = deseq2_values, edger_3_500_500_DOWNREGULATED = edger_values, noiseq_3_500_500_DOWNREGULATED = noiseq_values)


library(VennDiagram)

# Create the Venn diagram
venn.diagram(
  x = list(
    deseq2 = deseq2_values[!is.na(deseq2_values)],
    edger = edger_values[!is.na(edger_values)],
    noiseq = noiseq_values[!is.na(noiseq_values)]
  ),
  category.names = c("deseq2", "edger", "noiseq"),
  output = TRUE,
  filename = "Working Directory/Output/Images/venn_3_500_500_downregulated.png",
  output.type = "png",
  imagetype = "png",
  resolution = 300,
  category.col = c("red", "blue", "green"), #btw im colourblind so please feel free to update the colours
  fill = c("red", "blue", "green") 
  
)