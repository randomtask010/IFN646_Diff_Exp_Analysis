run_loop_noiseq_threshold <-function(SourceFileVariable, q_value) {
  
  # Step 1: Call Library (NOISeq)
  library(NOISeq)
  
  # Step 2: Load Dataset
  count_data <- read.table(paste0("RAW data/", SourceFileVariable, ".tsv"), header=TRUE, row.names=1)
  
  # Step 3: Create a conditions factor and a matrix of data
  samplesize <- ncol(count_data) / 2
  condition <- rep(c("condition1", "condition2"), each = samplesize)
  sample_names <- paste(condition, 1:(2*samplesize), sep="_") 
  colnames(count_data) <- sample_names
  sampleinfo <- data.frame(row.names=sample_names, condition=condition)
  
  # Creating a data frame for the sample conditions
  sampleinfo <- data.frame(row.names=sample_names, condition=condition)
  sampleinfo
  
  #Filtering data
  #Excluding features with low counts improves differential expression results, since noise in the data is reduced
  #counts per million (CPM) method is used to filtering
  count_data = filtered.data(count_data, factor = sampleinfo$condition, norm=FALSE, depth=NULL, method=1, cv.cutoff=100, cpm=1, p.adj="fdr")
  count_data
  
  # Step 4: Create NOISeq object
  mydata <- readData(data=count_data, factors=sampleinfo)
  mydata
  
  # Step 5: Differential expression analysis with NOISeq
  results <- noiseq(mydata, k = 0.5, norm = "tmm", replicates = "technical", factor = "condition")
  
  DEgenes <- degenes(results, q = q_value, M = NULL)  
  meta_data <- read.table(paste0("RAW data/", SourceFileVariable,"_meta.tsv"), header=TRUE, row.names=1)
  
  annotated_results <- merge(as.data.frame(DEgenes), meta_data, by="row.names", all.x=TRUE)
  annotated_results
  rownames(annotated_results) <- annotated_results$Row.names
  annotated_results$Row.names <- NULL
  
  # Downregulated in condition 2 compared to condition 1
  detected_down = degenes(results, q = q_value, M = "up")
  detected_down
  
  annotated_results_down <- merge(as.data.frame(detected_down), meta_data, by="row.names", all.x=TRUE)
  annotated_results_down
  rownames(annotated_results_down) <- annotated_results_down$Row.names
  annotated_results_down$Row.names <- NULL
  detected_down_NOISeq <- rownames(annotated_results_down)
  detected_down_NOISeq
  
  #from differentially expressed genes, actually downregulated according to metadata
  meta_down <- rownames(annotated_results[annotated_results$downregulation == 1,])
  
  common_down <- intersect(detected_down_NOISeq, meta_down)
  
  #detected_up represent the up regulated genes in condition 2 compared to condition 1.
  detected_up = degenes(results, q = q_value, M = "down")
  
  annotated_results_up <- merge(as.data.frame(detected_up), meta_data, by="row.names", all.x=TRUE)
  rownames(annotated_results_up) <- annotated_results_up$Row.names
  annotated_results_up$Row.names <- NULL
  detected_up_NOISeq <- rownames(annotated_results_up)
  
  
  #from differentially expressed genes, actually downregulated according to metadata
  meta_up <- rownames(annotated_results[annotated_results$upregulation == 1,])
  
  common_up <- intersect(detected_up_NOISeq, meta_up)
  
  # Step 10: Summarize outliers
  outliers_up <- setdiff(detected_up_NOISeq, meta_up)
  
  # Step 11: Accuracy and Precision Matrix
  true_positives <- length(common_up) + length(common_down)
  false_positives <- length(setdiff(detected_up_NOISeq, meta_up)) + length(setdiff(detected_down_NOISeq, meta_down))
  true_negatives <- nrow(meta_data) - (length(meta_up) + length(meta_down)) - false_positives
  false_negatives <- length(setdiff(meta_up, detected_up_NOISeq)) + length(setdiff(meta_down, detected_down_NOISeq))
  accuracy <- (true_positives + true_negatives) / nrow(meta_data)
  precision <- true_positives / (true_positives + false_positives)
  recall <- true_positives / (true_positives + false_negatives)
  f1_score <- 2 * ((precision * recall) / (precision + recall))
  FDR <- false_positives/(true_positives+false_positives) ## Precision = 1- FDR 
  
  # Step 12: Output Metrics to CSV for further analysis outside R
  metrics_df_1 <- data.frame(
    Threshold = q_value,
    True_Positives = true_positives,
    False_Positives = false_positives,
    True_Negatives = true_negatives,
    False_Negatives = false_negatives,
    Accuracy = accuracy,
    Precision = precision,
    Recall = recall,
    F1_Score = f1_score,
    FDR = FDR,
    Experiment = SourceFileVariable
  )

  # Initialize an empty data frame to hold all the results
  all_results <- data.frame()
  
  # Read the existing results, if any
  if (file.exists("Working Directory/Output/Threshold_noiseq.csv")) {
    all_results <- read.csv("Working Directory/Output/Threshold_noiseq.csv", header = TRUE)
  }
  
  # Append the current results to the all_results data frame
  all_results <- rbind(all_results, metrics_df_1)
  
  # Write the combined data frame to the CSV file
  write.csv(all_results, "Working Directory/Output/Threshold_noiseq.csv", row.names = FALSE)

}
