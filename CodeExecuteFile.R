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


  #parametrerised attempt - in progress
  source("Working Directory/Data_Analysis/Metrics_Dataframe.R")
  
  # down regulated stitch together for 3_500_500 down regulated outlier genes
  source("Working Directory/Data_Analysis/POC_3_500_500_DownRegulated.R")
  

  

#Analysis outputs and pretty pictures
  
  #pretty bar graph of metrics dataframe
  source("Working Directory/Data_Analysis/Accuracy_Metrics_Bar_Plot.R")
  