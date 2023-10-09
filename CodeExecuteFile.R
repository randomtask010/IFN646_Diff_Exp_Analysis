# Main Run Script

# Run Analysis

  # Run DESeq2
  
  source("Working Directory/Deseq2_Analysis/deseq2_3_500_500.R")
  source("Working Directory/Deseq2_Analysis/deseq2_3_750_250.R")
  source("Working Directory/Deseq2_Analysis/deseq2_3_1000_0.R")
  
  source("Working Directory/Deseq2_Analysis/deseq2_6_500_500.R")
  source("Working Directory/Deseq2_Analysis/deseq2_6_750_250.R")
  source("Working Directory/Deseq2_Analysis/deseq2_6_1000_0.R")
  
  source("Working Directory/Deseq2_Analysis/deseq2_9_500_500.R")
  source("Working Directory/Deseq2_Analysis/deseq2_9_750_250.R")
  source("Working Directory/Deseq2_Analysis/deseq2_9_1000_0.R")
  
  
  #Run Edge R Files
  source("Working Directory/EdgeR_Analysis/edgeR_3_500_500.R")
  source("Working Directory/EdgeR_Analysis/edgeR_3_750_250.R")
  source("Working Directory/EdgeR_Analysis/edgeR_3_1000_0.R")
  
  source("Working Directory/EdgeR_Analysis/edgeR_6_500_500.R")
  source("Working Directory/EdgeR_Analysis/edgeR_6_750_250.R")
  source("Working Directory/EdgeR_Analysis/edgeR_6_1000_0.R")
  
  source("Working Directory/EdgeR_Analysis/edgeR_9_500_500.R")
  source("Working Directory/EdgeR_Analysis/edgeR_9_750_250.R")
  source("Working Directory/EdgeR_Analysis/edgeR_9_1000_0.R")
  
  # Run NOISeq
  
  source("Working Directory/Noiseq_Analysis/noiseq_3_500_500.R")
  source("Working Directory/Noiseq_Analysis/noiseq_3_750_250.R")
  source("Working Directory/Noiseq_Analysis/noiseq_3_1000_0.R")
  
  source("Working Directory/Noiseq_Analysis/noiseq_6_500_500.R")
  source("Working Directory/Noiseq_Analysis/noiseq_6_750_250.R")
  source("Working Directory/Noiseq_Analysis/noiseq_6_1000_0.R")
  
  source("Working Directory/Noiseq_Analysis/noiseq_9_500_500.R")
  source("Working Directory/Noiseq_Analysis/noiseq_9_750_250.R")
  source("Working Directory/Noiseq_Analysis/noiseq_9_1000_0.R")


#Load in Metrics into a combined dataframe


  # POC for 3_500_500 - THIS WILL NOT STAY, but its a 1am POC that will need to be paramertised and made elegant - we will eventually move away from the manual csv comparisons for up/down regulated comparisons of overlaping mis-identified genes
  desq2_3_500_500_DOWN_df <- read.csv("Working Directory/Output/deseq2_3_500_500_outliers_downregulated.csv", stringsAsFactors = FALSE)
  edger_3_500_500_DOWN_df <- read.csv("Working Directory/Output/edger_3_500_500_outliers_downregulated.csv", stringsAsFactors = FALSE)
  noiseq_3_500_500_DOWN_df <- read.csv("Working Directory/Output/noiseq_3_500_500_outliers_downregulated.csv", stringsAsFactors = FALSE)
  deseq2_values <- desq2_3_500_500_DOWN_df[[1]]
  edger_values <- edger_3_500_500_DOWN_df[[1]]
  noiseq_values <- noiseq_3_500_500_DOWN_df[[1]]
  max_rows <- max(length(deseq2_values), length(edger_values), length(noiseq_values))
  deseq2_values <- c(deseq2_values, rep(NA, max_rows - length(deseq2_values)))
  edger_values <- c(edger_values, rep(NA, max_rows - length(edger_values)))
  noiseq_values <- c(noiseq_values, rep(NA, max_rows - length(noiseq_values)))
  downregulated_genes_3_500_500 <- data.frame(desq2_3_500_500_DOWNREGULATED = deseq2_values, edger_3_500_500_DOWNREGULATED = edger_values, noiseq_3_500_500_DOWNREGULATED = noiseq_values)

  #parametrerised attempt - in progress
  OutputDir <- "Working Directory/Output/"
  DownregulatedStr <- "_outliers_downregulated.csv"
  UpRegulatedStr <- "_outliers_upregulated.csv"
  groupings <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")
  toolset <- c("deseq2", "edger", "noiseq")
  
  
  

  

#Analysis outputs and pretty pictures
  
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
    filename = "Working Directory/Output/venn_3_500_500_downregulated.png",
    output.type = "png",
    imagetype = "png",
    resolution = 300,
    category.col = c("red", "blue", "green"), #btw im colourblind so please feel free to update the colours
    fill = c("red", "blue", "green") 
  
  )
  