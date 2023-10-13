# Main Run Script

# Run Analysis

  # Run DESeq2
  source("Working Directory/Deseq2_Analysis/deseq2_parameterised_file.R")
  #Variables
  Tool <- "deseq2_baseline"
  SourceFileVariable <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")
  PValue <- 0.05
  # loop
  for (sample in SourceFileVariable){
    run_loop_deseq2(Tool, sample, PValue)
  }
  
  
  #Run Edge R Files
  source("Working Directory/EdgeR_Analysis/edgeR_parameterised_file.R")
  #Variables
  Tool <- "edgeR_baseline"
  SourceFileVariable <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")
  PValue <- 0.05
  # loop
  for (sample in SourceFileVariable){
    run_loop_edgeR(Tool, sample, PValue)
  }
  
  # Run NOISeq
  source("Working Directory/Noiseq_Analysis/noiseq_parameterised_file.R")
  #Variables
  Tool <- "noiseq_baseline"
  SourceFileVariable <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")
  QValue <- 0.8
  # loop
  for (sample in SourceFileVariable){
    run_loop_noiseq(Tool, sample, QValue)
  }
  


# Metrics and Analysis


  #parametrerised attempt - in progress
  source("Working Directory/Data_Analysis/Metrics_Dataframe.R")
  
  # down regulated stitch together for 3_500_500 down regulated outlier genes
  source("Working Directory/Data_Analysis/POC_3_500_500_DownRegulated.R")
  
  ## Threshold Analysis ###
  # Run DeSeq2 Threshold analysis for all the files and save it into newly created Threshold_deseq.csv
  if (file.exists("Working Directory/Output/Threshold_deseq.csv")) {
    file.remove("Working Directory/Output/Threshold_deseq.csv")
  }
  source("Working Directory/Data_Analysis/Threshold_Analysis_Deseq_Parameterized.R")
  #Variables
  SourceFileVariable <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")
  PValue <- c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09)
  Tool <- "deseq2"
  # loop
  for (sample in SourceFileVariable){
    for(p in PValue )
    {
      run_loop_deseq_threshold(sample, p)
    }
  }
  
  # Run Edger Threshold analysis for all the files and save it into newly created Threshold_edger.csv
  if (file.exists("Working Directory/Output/Threshold_edger.csv")) {
    file.remove("Working Directory/Output/Threshold_edger.csv")
  }
  source("Working Directory/Data_Analysis/Threshold_Analysis_Edger_Parameterized.R")
  #Variables
  SourceFileVariable <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")
  PValue <- c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09)
  Tool <- "edgeR"
  # loop
  for (sample in SourceFileVariable){
    for(p in PValue )
    {
      run_loop_edgeR_threshold(sample, p)
    }
  }
  
  
  # Run NOISeq Threshold analysis for all the files and save it into newly created Threshold_noiseq.csv
  if (file.exists("Working Directory/Output/Threshold_noiseq.csv")) {
    file.remove("Working Directory/Output/Threshold_noiseq.csv")
  }
  source("Working Directory/Data_Analysis/Threshold_Analysis_Noiseq_Parameterized.R")
  #Variables
  SourceFileVariable <- c("3_500_500", "3_750_250", "3_1000_0", "6_500_500", "6_750_250", "6_1000_0", "9_500_500", "9_750_250", "9_1000_0")
  QValue <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
  Tool <- "noiseq"
  # loop
  for (sample in SourceFileVariable){
    for(q in QValue )
    {
      run_loop_noiseq_threshold(sample, q)
    }
  }
  
  
  
  # Outliers from Down Regulation
  source("Working Directory/Data_Analysis/Outliers_downregulation.R")
  PValue <- c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09)
  QValue <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
  # loop
      for(p in PValue ){
        for(q in QValue)
          run_loop_DownRegulation(p,q)
    }
  
  
  
  # Outliers from Up  Regulation
  source("Working Directory/Data_Analysis/Outliers_upregulation.R")
  PValue <- c(0.01,0.02,0.03,0.04,0.05,0.06,0.07,0.08,0.09)
  QValue <- c(0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9)
  # loop
  for(p in PValue ){
    for(q in QValue)
      run_loop_DownRegulation(p,q)
  }
  
# Analysis outputs and pretty pictures
  
  #pretty bar graph of metrics dataframe
  source("Working Directory/Data_Analysis/Confusion_Metrices_Bar_Plot.R")
  source("Working Directory/Data_Analysis/SampleSize_Impact_Analysis.R")
  source("Working Directory/Data_Analysis/GeneCount_Analysis_Samplesize3.R")
  source("Working Directory/Data_Analysis/GeneCount_Analysis_Samplesize6.R")
  source("Working Directory/Data_Analysis/GeneCount_Analysis_Samplesize9.R")
  source("Working Directory/Data_Analysis/Threshold_Analysis_Noiseq.R")
  source("Working Directory/Data_Analysis/Threshold_Analysis_Deseq.R")
  source("Working Directory/Data_Analysis/Threshold_Analysis_Edger.R")
  
  
  
  
  