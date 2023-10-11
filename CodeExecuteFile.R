# Main Run Script

# Run Analysis

  # Run DESeq2
  source("Working Directory/Deseq2_Analysis/deseq2_parameterised_file.R")
  #Variables
  Tool <- "deseq2"
  SourceFileVariable <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")
  PValue <- 0.05
  # loop
  for (sample in SourceFileVariable){
    run_loop_deseq2(Tool, sample, PValue)
  }
  
  
  #Run Edge R Files
  source("Working Directory/EdgeR_Analysis/edgeR_parameterised_file.R")
  #Variables
  Tool <- "edgeR"
  SourceFileVariable <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")
  PValue <- 0.05
  # loop
  for (sample in SourceFileVariable){
    run_loop_edgeR(Tool, sample, PValue)
  }
  
  # Run NOISeq
  source("Working Directory/Noiseq_Analysis/noiseq_parameterised_file.R")
  #Variables
  Tool <- "noiseq"
  SourceFileVariable <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")
  QValue <- 0.08
  # loop
  for (sample in SourceFileVariable){
    run_loop_noiseq(Tool, sample, QValue)
  }
  


#Load in Metrics into a combined dataframe


  #parametrerised attempt - in progress
  source("Working Directory/Data_Analysis/Metrics_Dataframe.R")
  
  # down regulated stitch together for 3_500_500 down regulated outlier genes
  source("Working Directory/Data_Analysis/POC_3_500_500_DownRegulated.R")
  
  #NoiSeq Data Merge
  source("Working Directory/Data_Analysis/Threshold_Analysis_Noiseq_3_1000_0.R")
  source("Working Directory/Data_Analysis/Threshold_Analysis_Noiseq_6_1000_0.R")
  source("Working Directory/Data_Analysis/Threshold_Analysis_Noiseq_9_1000_0.R")
  
  #DeSeq Data Merge
  source("Working Directory/Data_Analysis/Threshold_Analysis_Deseq_3_1000_0.R")
  source("Working Directory/Data_Analysis/Threshold_Analysis_Deseq_6_1000_0.R")
  source("Working Directory/Data_Analysis/Thershold_Analysis_Deseq_9_1000_0.R")
  
  #Edger Data Merge
  source("Working Directory/Data_Analysis/Threshold_Analysis_Edger_3_1000_0.R")
  source("Working Directory/Data_Analysis/Threshold_Analysis_Edger_6_1000_0.R")
  source("Working Directory/Data_Analysis/Threshold_Analysis_Edger_9_1000_0.R")
  

#Analysis outputs and pretty pictures
  
  #pretty bar graph of metrics dataframe
  source("Working Directory/Data_Analysis/Accuracy_Metrics_Bar_Plot.R")
  source("Working Directory/Data_Analysis/SampleSize_Impact_Analysis.R")
  source("Working Directory/Data_Analysis/GeneCount_Analysis_Samplesize3.R")
  source("Working Directory/Data_Analysis/GeneCount_Analysis_Samplesize6.R")
  source("Working Directory/Data_Analysis/GeneCount_Analysis_Samplesize9.R")
  source("Working Directory/Data_Analysis/Threshold_Analysis_Noiseq.R")
  source("Working Directory/Data_Analysis/Threshold_Analysis_Deseq.R")
  source("Working Directory/Data_Analysis/Threshold_Analysis_Edger.R")
  
  
  
  
  
  